/**
 * @file ForwardConfigurator.cpp
 * @brief Source of the implementation of the base forward configurator in xD space
 */

// Configuration includes
//
#include "QuICC/Debug/DebuggerMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/TransformConfigurators/ForwardConfigurator.hpp"

// Project includes
//
#include "QuICC/Arithmetics/SetNeg.hpp"
#include "QuICC/Arithmetics/None.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"

namespace QuICC {

namespace Transform {

   void ForwardConfigurator::nonlinearTerm(const TransformTree& tree, Physical::Kernel::SharedIPhysicalKernel spKernel, TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::NONLINEAR);

      // Get physical storage
      auto pNLComp = coord.ss().fwdPtr(static_cast<Dimensions::Transform::Id>(coord.ss().dimension()-1));
      coord.communicator().providePhysical(pNLComp);

      // Compute physical space kernel
      std::visit([&](auto&& pNL){spKernel->compute(*pNL, tree.comp<FieldComponents::Physical::Id>());}, pNLComp);

      // Transfer physical storage to next step
      coord.communicator().holdPhysical(pNLComp);

      ProfilerMacro_stop(Debug::Profiler::NONLINEAR);
   }

   void ForwardConfigurator::integrateND(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate ND", 4);

      ProfilerMacro_start(Debug::Profiler::FWDND);
      ProfilerMacro_start(Debug::Profiler::FWDNDIN);

      Dimensions::Transform::Id traND = static_cast<Dimensions::Transform::Id>(coord.ss().dimension()-1);

      // Get recover input data from hold
      auto pInVar = coord.ss().fwdPtr(traND);
      coord.communicator().receiveForward(traND, pInVar);

      // Get output storage
      auto pOutVar = coord.ss().bwdPtr(traND);
      auto pRecOutVar = coord.ss().bwdPtr(traND);
      if(edge.recoverOutId() >= 0)
      {
         coord.communicator().storage(traND).provideBwd(pOutVar);
         coord.communicator().storage(traND).recoverBwd(pRecOutVar, edge.recoverOutId());
      } else
      {
         coord.communicator().storage(traND).provideBwd(pOutVar);
      }

      ProfilerMacro_stop(Debug::Profiler::FWDNDIN);
      ProfilerMacro_start(Debug::Profiler::FWDNDTRA);

      // Compute integration transform of third dimension
      std::visit([&](auto&& pOut, auto&& pIn){coord.transformND().forward(pOut->rData(), pIn->data(), edge.opId());}, pOutVar, pInVar);

      ProfilerMacro_stop(Debug::Profiler::FWDNDTRA);
      ProfilerMacro_start(Debug::Profiler::FWDNDOUT);
      ProfilerMacro_start(Debug::Profiler::FWDNDOUTCOMM);

      // Hold temporary input storage
      if(edge.holdInput())
      {
         coord.communicator().storage(traND).holdFwd(pInVar);

      // Free temporary input storage
      } else
      {
         coord.communicator().storage(traND).freeFwd(pInVar);
      }

      ProfilerMacro_stop(Debug::Profiler::FWDNDOUTCOMM);

      // Combine recovered output with new calculation
      if(std::visit([](auto&& pRecOut){return (pRecOut != 0);},pRecOutVar))
      {
         ProfilerMacro_start(Debug::Profiler::FWDNDOUTWORK);

         std::visit([&](auto&& pRecOut, auto&& pOut){Datatypes::FieldTools::combine(*pRecOut, *pOut, edge.combinedArithId());}, pRecOutVar, pOutVar);

         ProfilerMacro_stop(Debug::Profiler::FWDNDOUTWORK);
         ProfilerMacro_start(Debug::Profiler::FWDNDOUTCOMM);

         if(edge.combinedOutId() >= 0)
         {
            coord.communicator().storage(traND).holdBwd(pRecOutVar, edge.combinedOutId());
         } else
         {
            coord.communicator().transferBackward(traND, pRecOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::FWDNDOUTCOMM);

      // Hold data for combination
      } else if(edge.combinedOutId() >= 0)
      {
         ProfilerMacro_start(Debug::Profiler::FWDNDOUTWORK);

         if(edge.combinedArithId() == Arithmetics::SetNeg::id())
         {
            std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::FWDNDOUTWORK);
         ProfilerMacro_start(Debug::Profiler::FWDNDOUTCOMM);

         coord.communicator().storage(traND).holdBwd(pOutVar, edge.combinedOutId());

         ProfilerMacro_stop(Debug::Profiler::FWDNDOUTCOMM);
      }

      // Transfer calculation
      if(edge.arithId() != Arithmetics::None::id())
      {
         ProfilerMacro_start(Debug::Profiler::FWDNDOUTWORK);

         if(edge.arithId() == Arithmetics::SetNeg::id())
         {
            std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::FWDNDOUTWORK);
         ProfilerMacro_start(Debug::Profiler::FWDNDOUTCOMM);

         coord.communicator().transferBackward(traND, pOutVar);

         ProfilerMacro_stop(Debug::Profiler::FWDNDOUTCOMM);

      } else if(edge.combinedOutId() < 0 || std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
      {
         ProfilerMacro_start(Debug::Profiler::FWDNDOUTCOMM);

         coord.communicator().storage(traND).freeBwd(pOutVar);

         ProfilerMacro_stop(Debug::Profiler::FWDNDOUTCOMM);
      }

      ProfilerMacro_stop(Debug::Profiler::FWDNDOUT);
      ProfilerMacro_stop(Debug::Profiler::FWDND);
   }

   void ForwardConfigurator::integrate2D(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate 2D", 4);

      ProfilerMacro_start(Debug::Profiler::FWD2D);
      ProfilerMacro_start(Debug::Profiler::FWD2DIN);

      auto pInVar = coord.ss().fwdPtr(Dimensions::Transform::TRA2D);
      if(edge.recoverInput())
      {
         coord.communicator().storage(Dimensions::Transform::TRA2D).recoverFwd(pInVar);

      // Get the transfered input data
      } else
      {
         coord.communicator().receiveForward(Dimensions::Transform::TRA2D, pInVar);
      }

      // Get output storage
      auto pOutVar = coord.ss().bwdPtr(Dimensions::Transform::TRA2D);
      auto pRecOutVar = coord.ss().bwdPtr(Dimensions::Transform::TRA2D);
      if(edge.recoverOutId() >= 0)
      {
         coord.communicator().storage(Dimensions::Transform::TRA2D).provideBwd(pOutVar);
         coord.communicator().storage(Dimensions::Transform::TRA2D).recoverBwd(pRecOutVar, edge.recoverOutId());
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA2D).provideBwd(pOutVar);
      }

      ProfilerMacro_stop(Debug::Profiler::FWD2DIN);
      ProfilerMacro_start(Debug::Profiler::FWD2DTRA);

      // Compute integration transform of second dimension
      std::visit([&](auto&& pOut, auto&& pIn){coord.transform2D().forward(pOut->rData(), pIn->data(), edge.opId());}, pOutVar, pInVar);

      ProfilerMacro_stop(Debug::Profiler::FWD2DTRA);
      ProfilerMacro_start(Debug::Profiler::FWD2DOUT);
      ProfilerMacro_start(Debug::Profiler::FWD2DOUTCOMM);

      // Hold temporary storage
      if(edge.holdInput())
      {
         coord.communicator().storage(Dimensions::Transform::TRA2D).holdFwd(pInVar);

      // Free temporary storage
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA2D).freeFwd(pInVar);
      }

      ProfilerMacro_stop(Debug::Profiler::FWD2DOUTCOMM);

      // Combine recovered output with new calculation
      if(std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
      {
         ProfilerMacro_start(Debug::Profiler::FWD2DOUTWORK);

         std::visit([&](auto&& pRecOut, auto&& pOut){Datatypes::FieldTools::combine(*pRecOut, *pOut, edge.combinedArithId());}, pRecOutVar, pOutVar);

         ProfilerMacro_stop(Debug::Profiler::FWD2DOUTWORK);
         ProfilerMacro_start(Debug::Profiler::FWD2DOUTCOMM);

         if(edge.combinedOutId() >= 0)
         {
            coord.communicator().storage(Dimensions::Transform::TRA2D).holdBwd(pRecOutVar, edge.combinedOutId());
         } else
         {
            coord.communicator().transferBackward(Dimensions::Transform::TRA2D, pRecOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::FWD2DOUTCOMM);

      // Hold data for combination
      } else if(edge.combinedOutId() >= 0)
      {
         ProfilerMacro_start(Debug::Profiler::FWD2DOUTWORK);

         if(edge.combinedArithId() == Arithmetics::SetNeg::id())
         {
            std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::FWD2DOUTWORK);
         ProfilerMacro_start(Debug::Profiler::FWD2DOUTCOMM);

         coord.communicator().storage(Dimensions::Transform::TRA2D).holdBwd(pOutVar, edge.combinedOutId());

         ProfilerMacro_stop(Debug::Profiler::FWD2DOUTCOMM);
      }

      // Transfer calculation
      if(edge.arithId() != Arithmetics::None::id())
      {
         ProfilerMacro_start(Debug::Profiler::FWD2DOUTWORK);

         if(edge.arithId() == Arithmetics::SetNeg::id())
         {
            std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::FWD2DOUTWORK);
         ProfilerMacro_start(Debug::Profiler::FWD2DOUTCOMM);

         coord.communicator().transferBackward(Dimensions::Transform::TRA2D, pOutVar);

         ProfilerMacro_stop(Debug::Profiler::FWD2DOUTCOMM);

      } else if(edge.combinedOutId() < 0 || std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
      {
         ProfilerMacro_start(Debug::Profiler::FWD2DOUTCOMM);

         coord.communicator().storage(Dimensions::Transform::TRA2D).freeBwd(pOutVar);

         ProfilerMacro_stop(Debug::Profiler::FWD2DOUTCOMM);
      }

      ProfilerMacro_stop(Debug::Profiler::FWD2DOUT);
      ProfilerMacro_stop(Debug::Profiler::FWD2D);
   }

   void ForwardConfigurator::integrate1D(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate 1D", 4);

      ProfilerMacro_start(Debug::Profiler::FWD1D);
      ProfilerMacro_start(Debug::Profiler::FWD1DIN);

      auto pInVar = coord.ss().fwdPtr(Dimensions::Transform::TRA1D);

      // Recover hold input data
      if(edge.recoverInput())
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).recoverFwd(pInVar);

      // Get the transfered input data
      } else
      {
         coord.communicator().receiveForward(Dimensions::Transform::TRA1D, pInVar);
      }

      // Get temporary storage
      auto pOutVar = coord.ss().bwdPtr(Dimensions::Transform::TRA1D);
      coord.communicator().storage(Dimensions::Transform::TRA1D).provideBwd(pOutVar);

      ProfilerMacro_stop(Debug::Profiler::FWD1DIN);
      ProfilerMacro_start(Debug::Profiler::FWD1DTRA);

      // Compute integration transform of first dimension
      std::visit([&](auto&& pOut, auto&& pIn){coord.transform1D().forward(pOut->rData(), pIn->data(), edge.opId());}, pOutVar, pInVar);

      ProfilerMacro_stop(Debug::Profiler::FWD1DTRA);
      ProfilerMacro_start(Debug::Profiler::FWD1DOUT);
      ProfilerMacro_start(Debug::Profiler::FWD1DOUTCOMM);

      // Hold temporary storage
      if(edge.holdInput())
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).holdFwd(pInVar);

      // Free temporary storage
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).freeFwd(pInVar);
      }

      // Hold temporary storage
      coord.communicator().storage(Dimensions::Transform::TRA1D).holdBwd(pOutVar);

      ProfilerMacro_stop(Debug::Profiler::FWD1DOUTCOMM);
      ProfilerMacro_stop(Debug::Profiler::FWD1DOUT);
      ProfilerMacro_stop(Debug::Profiler::FWD1D);
   }

   void ForwardConfigurator::integrate1ND(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate 1D", 4);

      ProfilerMacro_start(Debug::Profiler::FWD1D);
      ProfilerMacro_start(Debug::Profiler::FWDIN);

      // Get recover input data from hold
      auto pInVar = coord.ss().fwdPtr(Dimensions::Transform::TRA1D);
      coord.communicator().receiveForward(Dimensions::Transform::TRA1D, pInVar);

      // Get output storage
      auto pOutVar = coord.ss().bwdPtr(Dimensions::Transform::TRA1D);
      auto pRecOutVar = coord.ss().bwdPtr(Dimensions::Transform::TRA1D);
      if(edge.recoverOutId() >= 0)
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).provideBwd(pOutVar);
         coord.communicator().storage(Dimensions::Transform::TRA1D).recoverBwd(pRecOutVar, edge.recoverOutId());
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).provideBwd(pOutVar);
      }

      ProfilerMacro_stop(Debug::Profiler::FWD1DIN);
      ProfilerMacro_start(Debug::Profiler::FWD1DTRA);

      // Compute integration transform of third dimension
      std::visit([&](auto&& pOut, auto&& pIn){coord.transform1D().forward(pOut->rData(), pIn->data(), edge.opId());}, pOutVar, pInVar);

      ProfilerMacro_stop(Debug::Profiler::FWD1DTRA);
      ProfilerMacro_start(Debug::Profiler::FWD1DOUT);
      ProfilerMacro_start(Debug::Profiler::FWD1DOUTCOMM);

      // Hold temporary storage
      if(edge.holdInput())
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).holdFwd(pInVar);

      // Free temporary storage
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).freeFwd(pInVar);
      }

      // Hold temporary storage
      coord.communicator().storage(Dimensions::Transform::TRA1D).holdBwd(pOutVar);

      ProfilerMacro_stop(Debug::Profiler::FWD1DOUTCOMM);
      ProfilerMacro_stop(Debug::Profiler::FWD1DOUT);
      ProfilerMacro_stop(Debug::Profiler::FWD1D);
   }

   void ForwardConfigurator::updateEquation(const TransformTreeEdge& edge, Framework::Selector::VariantSharedScalarVariable& rScalar, TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::FWDOUT);

      // Recover temporary storage
      auto pInVar = coord.ss().bwdPtr(Dimensions::Transform::TRA1D);
      coord.communicator().storage(Dimensions::Transform::TRA1D).recoverBwd(pInVar);

      // Push variable into data pool queue
      std::visit([&](auto&& p){coord.communicator().storage(Dimensions::Transform::TRA1D).holdBwd(p->rDom(0).rPerturbation());}, rScalar);

      // Compute linear term component
      coord.communicator().updateSpectral(pInVar, edge.arithId());

      // Hold temporary storage
      if(edge.combinedOutId() >= 0)
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).holdBwd(pInVar, edge.combinedOutId());

      // Free the temporary storage
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).freeBwd(pInVar);
      }

      ProfilerMacro_stop(Debug::Profiler::FWDOUT);
   }

   void ForwardConfigurator::updateEquation(const TransformTreeEdge& edge, Framework::Selector::VariantSharedVectorVariable& rVector, TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::FWDOUT);

      // Recover temporary storage
      auto pInVar = coord.ss().bwdPtr(Dimensions::Transform::TRA1D);
      coord.communicator().storage(Dimensions::Transform::TRA1D).recoverBwd(pInVar);

      // Push variable into data pool queue
      std::visit([&](auto&& p){coord.communicator().storage(Dimensions::Transform::TRA1D).holdBwd(p->rDom(0).rPerturbation().rComp(edge.outId<FieldComponents::Spectral::Id>()));}, rVector);

      // Compute linear term component
      coord.communicator().updateSpectral(pInVar, edge.arithId());

      // Hold temporary storage
      if(edge.combinedOutId() >= 0)
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).holdBwd(pInVar, edge.combinedOutId());

      // Free the temporary storage
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).freeBwd(pInVar);
      }

      ProfilerMacro_stop(Debug::Profiler::FWDOUT);
   }

}
}
