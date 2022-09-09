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
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   void ForwardConfigurator::nonlinearTerm(const TransformTree& tree, Physical::Kernel::SharedIPhysicalKernel spKernel, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("nonlinearTerm", 4);

      Profiler::RegionFixture<2> fix("Fwd-nonlinearTerm");

      // Get physical storage
      auto pNLComp = coord.ss().fwdPtr(static_cast<Dimensions::Transform::Id>(coord.ss().dimension()-1));
      coord.communicator().providePhysical(pNLComp);

      // Compute physical space kernel
      std::visit([&](auto&& pNL){spKernel->compute(*pNL, tree.comp<FieldComponents::Physical::Id>());}, pNLComp);

      // Transfer physical storage to next step
      coord.communicator().holdPhysical(pNLComp);
   }

   void ForwardConfigurator::integrateND(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate ND", 4);

      const std::string profRegion = "Fwd-integrateND";
      const auto traId = static_cast<Dimensions::Transform::Id>(coord.ss().dimension()-1);
      const bool processOutput = true;
      ForwardConfigurator::genericIntegrate(edge, coord, traId, processOutput, profRegion);
   }

   void ForwardConfigurator::integrate2D(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate 2D", 4);

      const std::string profRegion = "Fwd-integrate2D";
      const auto traId = Dimensions::Transform::TRA2D;
      const bool processOutput = true;
      ForwardConfigurator::genericIntegrate(edge, coord, traId, processOutput, profRegion);
   }

   void ForwardConfigurator::integrate1D(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate 1D", 4);

      const std::string profRegion = "Fwd-integrate1D";
      const auto traId = Dimensions::Transform::TRA1D;
      const bool processOutput = false;
      ForwardConfigurator::genericIntegrate(edge, coord, traId, processOutput, profRegion);
   }

   void ForwardConfigurator::genericIntegrate(const TransformTreeEdge& edge, TransformCoordinatorType& coord, const Dimensions::Transform::Id traId, const bool processOutput, const std::string profRegion)
   {
      Profiler::RegionFixture<2> fix(profRegion);
      Profiler::RegionStart<3> (profRegion + "-pre");

      auto pInVar = coord.ss().fwdPtr(traId);

      // Recover hold input data
      if(edge.recoverInput())
      {
         coord.communicator().storage(traId).recoverFwd(pInVar);
      }
      // Get the transfered input data
      else
      {
         coord.communicator().receiveForward(traId, pInVar);
      }

      // Get output storage
      auto pOutVar = coord.ss().bwdPtr(traId);
      coord.communicator().storage(traId).provideBwd(pOutVar);

      Profiler::RegionStop<3> (profRegion + "-pre");
      Profiler::RegionStart<3> (profRegion + "-transform");

      // Compute integration transform
      std::visit(
            [&](auto&& pOut, auto&& pIn)
            {
               coord.transform(traId).forward(pOut->rData(), pIn->data(), edge.opId());
            },
            pOutVar, pInVar);

      Profiler::RegionStop<3> (profRegion + "-transform");
      Profiler::RegionStart<3> (profRegion + "-post");
      Profiler::RegionStart<4> (profRegion + "-post" + "-comm_a");

      // Hold temporary input storage
      if(edge.holdInput())
      {
         coord.communicator().storage(traId).holdFwd(pInVar);
      }
      // Free temporary input storage
      else
      {
         coord.communicator().storage(traId).freeFwd(pInVar);
      }

      Profiler::RegionStop<4> (profRegion + "-post" + "-comm_a");

      if(processOutput)
      {
         Profiler::RegionStart<4> (profRegion + "-post" + "-work_b");

         // Combine recovered output with new calculation
         auto pRecOutVar = coord.ss().bwdPtr(traId);
         if(edge.recoverOutId() >= 0)
         {
            coord.communicator().storage(traId).recoverBwd(pRecOutVar, edge.recoverOutId());
         }

         Profiler::RegionStop<4> (profRegion + "-post" + "-work_b");

         if(std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
         {
            Profiler::RegionStart<4> (profRegion + "-post" + "-work_c");

            std::visit([&](auto&& pRecOut, auto&& pOut){Datatypes::FieldTools::combine(*pRecOut, *pOut, edge.combinedArithId());}, pRecOutVar, pOutVar);

            Profiler::RegionStop<4> (profRegion + "-post" + "-work_c");
            Profiler::RegionStart<4> (profRegion + "-post" + "-comm_c");

            if(edge.combinedOutId() >= 0)
            {
               coord.communicator().storage(traId).holdBwd(pRecOutVar, edge.combinedOutId());
            }
            else
            {
               coord.communicator().transferBackward(traId, pRecOutVar);
            }

            Profiler::RegionStop<4> (profRegion + "-post" + "-comm_c");
         }
         // Hold data for combination
         else if(edge.combinedOutId() >= 0)
         {
            Profiler::RegionStart<4> (profRegion + "-post" + "-work_d");

            if(edge.combinedArithId() == Arithmetics::SetNeg::id())
            {
               std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
            }

            Profiler::RegionStop<4> (profRegion + "-post" + "-work_d");
            Profiler::RegionStart<4> (profRegion + "-post" + "-comm_d");

            coord.communicator().storage(traId).holdBwd(pOutVar, edge.combinedOutId());

            Profiler::RegionStop<4> (profRegion + "-post" + "-comm_d");
         }

         // Transfer calculation
         if(edge.arithId() != Arithmetics::None::id())
         {
            Profiler::RegionStart<4> (profRegion + "-post" + "-work_e");

            if(edge.arithId() == Arithmetics::SetNeg::id())
            {
               std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
            }

            Profiler::RegionStop<4> (profRegion + "-post" + "-work_e");
            Profiler::RegionStart<4> (profRegion + "-post" + "-comm_e");

            coord.communicator().transferBackward(traId, pOutVar);

            Profiler::RegionStop<4> (profRegion + "-post" + "-comm_e");

         }
         else if(edge.combinedOutId() < 0 || std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
         {
            Profiler::RegionStart<4> (profRegion + "-post" + "-comm_f");

            coord.communicator().storage(traId).freeBwd(pOutVar);

            Profiler::RegionStop<4> (profRegion + "-post" + "-comm_f");
         }
      }
      else
      {
         Profiler::RegionStart<4> (profRegion + "-post" + "-comm_g");

         // Transfer data
         coord.communicator().transferBackward(traId, pOutVar);

         Profiler::RegionStop<4> (profRegion + "-post" + "-comm_g");
      }

      Profiler::RegionStop<3> (profRegion + "-post");
   }

   void ForwardConfigurator::integrate1ND(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate 1D", 4);

      const std::string profRegion = "Fwd-integrate1ND";
      throw std::logic_error("Foward 1ND setup is not ported");
#if 0
      Profiler::RegionFixture<2> fix(profRegion);
      Profiler::RegionStart<3> (profRegion + "-pre");

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

      Profiler::RegionStop<3> (profRegion + "-pre");
      Profiler::RegionStart<3> (profRegion + "-transform");

      // Compute integration transform of third dimension
      std::visit([&](auto&& pOut, auto&& pIn){coord.transform1D().forward(pOut->rData(), pIn->data(), edge.opId());}, pOutVar, pInVar);

      Profiler::RegionStop<3> (profRegion + "-transform");
      Profiler::RegionStart<3> (profRegion + "-post");
      Profiler::RegionStart<4> (profRegion + "-post" + "-comm_a");

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

      Profiler::RegionStop<4> (profRegion + "-post" + "-comm_a");
      Profiler::RegionStop<3> (profRegion + "-post");
#endif
   }

   void ForwardConfigurator::updateEquation(const TransformTreeEdge& edge, Framework::Selector::VariantSharedScalarVariable& rScalar, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Fwd-updateEquation");

      // Debugger message
      DebuggerMacro_msg("updateEquation (scalar)", 4);

      const auto transId = Dimensions::Transform::TRA1D;

      // Recover temporary storage
      auto pInVar = coord.ss().bwdPtr(transId);
      coord.communicator().receiveBackward(transId, pInVar);

      // Push variable into data pool queue
      std::visit([&](auto&& p){coord.communicator().storage(transId).holdBwd(p->rDom(0).rPerturbation());}, rScalar);

      // Compute linear term component
      coord.communicator().updateSpectral(pInVar, edge.arithId());

      // Hold temporary storage
      if(edge.combinedOutId() >= 0)
      {
         coord.communicator().storage(transId).holdBwd(pInVar, edge.combinedOutId());

      // Free the temporary storage
      } else
      {
         coord.communicator().storage(transId).freeBwd(pInVar);
      }
   }

   void ForwardConfigurator::updateEquation(const TransformTreeEdge& edge, Framework::Selector::VariantSharedVectorVariable& rVector, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Fwd-updateEquation");

      // Debugger message
      DebuggerMacro_msg("updateEquation (vector)", 4);

      const auto transId = Dimensions::Transform::TRA1D;

      // Recover temporary storage
      auto pInVar = coord.ss().bwdPtr(transId);
      coord.communicator().receiveBackward(transId, pInVar);

      // Push variable into data pool queue
      std::visit([&](auto&& p){coord.communicator().storage(transId).holdBwd(p->rDom(0).rPerturbation().rComp(edge.outId<FieldComponents::Spectral::Id>()));}, rVector);

      // Compute linear term component
      coord.communicator().updateSpectral(pInVar, edge.arithId());

      // Hold temporary storage
      if(edge.combinedOutId() >= 0)
      {
         coord.communicator().storage(transId).holdBwd(pInVar, edge.combinedOutId());

      // Free the temporary storage
      } else
      {
         coord.communicator().storage(transId).freeBwd(pInVar);
      }
   }

}
}
