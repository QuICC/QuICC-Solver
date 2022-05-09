/**
 * @file BackwardConfigurator.cpp
 * @brief Source of the implementation of the base backward configurator
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
#include "QuICC/TransformConfigurators/BackwardConfigurator.hpp"

// Project includes
//
#include "QuICC/Arithmetics/SetNeg.hpp"
#include "QuICC/Arithmetics/None.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   void BackwardConfigurator::prepareSpectral(const TransformTree& tree, Framework::Selector::VariantSharedScalarVariable& rScalar, TransformCoordinatorType& coord)
   {
      // Safety assert
      assert(tree.comp<FieldComponents::Spectral::Id>() == FieldComponents::Spectral::SCALAR);

      ProfilerMacro_start(Debug::Profiler::BWDIN);

      // Put scalar into temporary hold storage
      std::visit([&](auto&& p){coord.communicator().dealiasSpectral(p->rDom(0).rTotal());}, rScalar);

      ProfilerMacro_stop(Debug::Profiler::BWDIN);
   }

   void BackwardConfigurator::prepareSpectral(const TransformTree& tree, Framework::Selector::VariantSharedVectorVariable& rVector, TransformCoordinatorType& coord)
   {
      // Safety assert
      assert(tree.comp<FieldComponents::Spectral::Id>() != FieldComponents::Spectral::SCALAR);

      ProfilerMacro_start(Debug::Profiler::BWDIN);

      // Put scalar into temporary hold storage
      std::visit([&](auto&& p){coord.communicator().dealiasSpectral(p->rDom(0).rTotal().rComp(tree.comp<FieldComponents::Spectral::Id>()));}, rVector);

      ProfilerMacro_stop(Debug::Profiler::BWDIN);
   }

   void BackwardConfigurator::preparePhysical(const TransformTree&, const TransformTreeEdge& edge, Framework::Selector::VariantSharedScalarVariable& rScalar, TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::BWDOUT);

      // Put scalar into temporary hold storage
      if(edge.fieldId() == FieldType::SCALAR)
      {
         std::visit([&](auto&& p){coord.communicator().holdPhysical(p->rDom(0).rPhys());}, rScalar);

      // Put gradient component into temporary hold storage
      } else if(edge.fieldId() == FieldType::GRADIENT)
      {
         std::visit([&](auto&& p){coord.communicator().holdPhysical(p->rDom(0).rGrad().rComp(edge.outId<FieldComponents::Physical::Id>()));}, rScalar);

      // Put 2nd order gradient component into temporary hold storage
      } else if(edge.fieldId() == FieldType::GRADIENT2)
      {
         std::visit([&](auto&& p){coord.communicator().holdPhysical(p->rDom(0).rGrad2().rComp(edge.outId<FieldComponents::Physical::Id>(0),edge.outId<FieldComponents::Physical::Id>(1)));}, rScalar);
      }

      ProfilerMacro_stop(Debug::Profiler::BWDOUT);
   }

   void BackwardConfigurator::preparePhysical(const TransformTree& tree, const TransformTreeEdge& edge, Framework::Selector::VariantSharedVectorVariable& rVector, TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::BWDOUT);

      // Put vector component into temporary hold storage
      if(edge.fieldId() == FieldType::VECTOR)
      {
         std::visit([&](auto&& p){coord.communicator().holdPhysical(p->rDom(0).rPhys().rComp(edge.outId<FieldComponents::Physical::Id>()));}, rVector);

      // Put vector gradient component into temporary hold storage
      } else if(edge.fieldId() == FieldType::GRADIENT)
      {
         std::visit([&](auto&& p){coord.communicator().holdPhysical(p->rDom(0).rGrad(tree.comp<FieldComponents::Spectral::Id>()).rComp(edge.outId<FieldComponents::Physical::Id>()));}, rVector);

      // Put curl component into temporary hold storage
      } else if(edge.fieldId() == FieldType::CURL)
      {
         std::visit([&](auto&& p){coord.communicator().holdPhysical(p->rDom(0).rCurl().rComp(edge.outId<FieldComponents::Physical::Id>()));}, rVector);
      }

      ProfilerMacro_stop(Debug::Profiler::BWDOUT);
   }

   void BackwardConfigurator::project1D(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("BwdProject1D");

      // Debugger message
      DebuggerMacro_msg("Project 1D", 4);

      ProfilerMacro_start(Debug::Profiler::BWD1D);
      ProfilerMacro_start(Debug::Profiler::BWD1DIN);

      // Get the input data from hold
      auto pInVar = coord.ss().bwdPtr(Dimensions::Transform::TRA1D);
      coord.communicator().storage(Dimensions::Transform::TRA1D).recoverBwd(pInVar);

      // Get output storage
      auto pOutVar = coord.ss().fwdPtr(Dimensions::Transform::TRA1D);
      auto pRecOutVar = coord.ss().fwdPtr(Dimensions::Transform::TRA1D);
      if(edge.recoverOutId() >= 0)
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).provideFwd(pOutVar);
         coord.communicator().storage(Dimensions::Transform::TRA1D).recoverFwd(pRecOutVar, edge.recoverOutId());
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).provideFwd(pOutVar);
      }

      ProfilerMacro_stop(Debug::Profiler::BWD1DIN);
      ProfilerMacro_start(Debug::Profiler::BWD1DTRA);

      // Compute projection transform for first dimension
      std::visit([&](auto&& pOut, auto&& pIn){coord.transform1D().backward(pOut->rData(), pIn->data(), edge.opId());}, pOutVar, pInVar);

      ProfilerMacro_stop(Debug::Profiler::BWD1DTRA);
      ProfilerMacro_start(Debug::Profiler::BWD1DOUT);
      ProfilerMacro_start(Debug::Profiler::BWD1DOUTCOMM);

      // Hold spectral input
      if(edge.holdInput())
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).holdBwd(pInVar);

      // Free spectral input
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).freeBwd(pInVar);
      }

      ProfilerMacro_stop(Debug::Profiler::BWD1DOUTCOMM);

      // Combine recovered output with new calculation
      if(std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
      {
         ProfilerMacro_start(Debug::Profiler::BWD1DOUTWORK);

         std::visit([&](auto&& pRecOut, auto&& pOut){Datatypes::FieldTools::combine(*pRecOut, *pOut, edge.combinedArithId());}, pRecOutVar, pOutVar);

         ProfilerMacro_stop(Debug::Profiler::BWD1DOUTWORK);
         ProfilerMacro_start(Debug::Profiler::BWD1DOUTCOMM);

         if(edge.combinedOutId() >= 0)
         {
            coord.communicator().storage(Dimensions::Transform::TRA1D).holdFwd(pRecOutVar, edge.combinedOutId());
         } else
         {
            coord.communicator().transferForward(Dimensions::Transform::TRA1D, pRecOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::BWD1DOUTCOMM);

      // Hold data for combination
      } else if(edge.combinedOutId() >= 0)
      {
         ProfilerMacro_start(Debug::Profiler::BWD1DOUTWORK);

         if(edge.combinedArithId() == Arithmetics::SetNeg::id())
         {
            std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::BWD1DOUTWORK);
         ProfilerMacro_start(Debug::Profiler::BWD1DOUTCOMM);

         coord.communicator().storage(Dimensions::Transform::TRA1D).holdFwd(pOutVar, edge.combinedOutId());

         ProfilerMacro_stop(Debug::Profiler::BWD1DOUTCOMM);
      }

      // Transfer calculation
      if(edge.arithId() != Arithmetics::None::id())
      {
         ProfilerMacro_start(Debug::Profiler::BWD1DOUTWORK);

         if(edge.arithId() == Arithmetics::SetNeg::id())
         {
            std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::BWD1DOUTWORK);
         ProfilerMacro_start(Debug::Profiler::BWD1DOUTCOMM);

         coord.communicator().transferForward(Dimensions::Transform::TRA1D, pOutVar);

         ProfilerMacro_stop(Debug::Profiler::BWD1DOUTCOMM);

      } else if(edge.combinedOutId() < 0 || std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
      {
         ProfilerMacro_start(Debug::Profiler::BWD1DOUTCOMM);

         coord.communicator().storage(Dimensions::Transform::TRA1D).freeFwd(pOutVar);

         ProfilerMacro_stop(Debug::Profiler::BWD1DOUTCOMM);
      }

      ProfilerMacro_stop(Debug::Profiler::BWD1DOUT);
      ProfilerMacro_stop(Debug::Profiler::BWD1D);
   }

   void BackwardConfigurator::project2D(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("BwdProject2D");

      // Debugger message
      DebuggerMacro_msg("Project 2D", 4);

      ProfilerMacro_start(Debug::Profiler::BWD2D);
      ProfilerMacro_start(Debug::Profiler::BWD2DIN);

      auto pInVar = coord.ss().bwdPtr(Dimensions::Transform::TRA2D);

      // Get the input data from hold
      if(edge.recoverInput())
      {
         coord.communicator().storage(Dimensions::Transform::TRA2D).recoverBwd(pInVar);

      // Get the transfered input data
      } else
      {
         coord.communicator().receiveBackward(Dimensions::Transform::TRA2D, pInVar);
      }

      // Get output storage
      auto pOutVar = coord.ss().fwdPtr(Dimensions::Transform::TRA2D);
      auto pRecOutVar = coord.ss().fwdPtr(Dimensions::Transform::TRA2D);
      if(edge.recoverOutId() >= 0)
      {
         coord.communicator().storage(Dimensions::Transform::TRA2D).provideFwd(pOutVar);
         coord.communicator().storage(Dimensions::Transform::TRA2D).recoverFwd(pRecOutVar, edge.recoverOutId());
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA2D).provideFwd(pOutVar);
      }

      ProfilerMacro_stop(Debug::Profiler::BWD2DIN);
      ProfilerMacro_start(Debug::Profiler::BWD2DTRA);

      // Compute projection transform for second dimension
      std::visit([&](auto&& pOut, auto&& pIn){coord.transform2D().backward(pOut->rData(), pIn->data(), edge.opId());}, pOutVar, pInVar);

      ProfilerMacro_stop(Debug::Profiler::BWD2DTRA);
      ProfilerMacro_start(Debug::Profiler::BWD2DOUT);
      ProfilerMacro_start(Debug::Profiler::BWD2DOUTCOMM);

      // Hold temporary storage
      if(edge.holdInput())
      {
         coord.communicator().storage(Dimensions::Transform::TRA2D).holdBwd(pInVar);

      // Free temporary input storage
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA2D).freeBwd(pInVar);
      }

      ProfilerMacro_stop(Debug::Profiler::BWD2DOUTCOMM);

      // Combine recovered output with new calculation
      if(std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
      {
         ProfilerMacro_start(Debug::Profiler::BWD2DOUTWORK);

         std::visit([&](auto&& pRecOut, auto&& pOut){Datatypes::FieldTools::combine(*pRecOut, *pOut, edge.combinedArithId());}, pRecOutVar, pOutVar);

         ProfilerMacro_stop(Debug::Profiler::BWD2DOUTWORK);
         ProfilerMacro_start(Debug::Profiler::BWD2DOUTCOMM);

         if(edge.combinedOutId() >= 0)
         {
            coord.communicator().storage(Dimensions::Transform::TRA2D).holdFwd(pRecOutVar, edge.combinedOutId());
         } else
         {
            coord.communicator().transferForward(Dimensions::Transform::TRA2D, pRecOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::BWD2DOUTCOMM);

      // Hold data for combination
      } else if(edge.combinedOutId() >= 0)
      {
         ProfilerMacro_start(Debug::Profiler::BWD2DOUTWORK);

         if(edge.combinedArithId() == Arithmetics::SetNeg::id())
         {
            std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::BWD2DOUTWORK);
         ProfilerMacro_start(Debug::Profiler::BWD2DOUTCOMM);

         coord.communicator().storage(Dimensions::Transform::TRA2D).holdFwd(pOutVar, edge.combinedOutId());

         ProfilerMacro_stop(Debug::Profiler::BWD2DOUTCOMM);
      }

      // Transfer calculation
      if(edge.arithId() != Arithmetics::None::id())
      {
         ProfilerMacro_start(Debug::Profiler::BWD2DOUTWORK);

         if(edge.arithId() == Arithmetics::SetNeg::id())
         {
            std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::BWD2DOUTWORK);
         ProfilerMacro_start(Debug::Profiler::BWD2DOUTCOMM);

         coord.communicator().transferForward(Dimensions::Transform::TRA2D, pOutVar);

         ProfilerMacro_stop(Debug::Profiler::BWD2DOUTCOMM);

      } else if(edge.combinedOutId() < 0 || std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
      {
         ProfilerMacro_start(Debug::Profiler::BWD2DOUTCOMM);

         coord.communicator().storage(Dimensions::Transform::TRA2D).freeFwd(pOutVar);

         ProfilerMacro_stop(Debug::Profiler::BWD2DOUTCOMM);
      }

      ProfilerMacro_stop(Debug::Profiler::BWD2DOUT);
      ProfilerMacro_stop(Debug::Profiler::BWD2D);
   }

   void BackwardConfigurator::projectND(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("BwdProjectND");

      // Debugger message
      DebuggerMacro_msg("Project ND", 4);

      ProfilerMacro_start(Debug::Profiler::BWDND);
      ProfilerMacro_start(Debug::Profiler::BWDNDIN);

      Dimensions::Transform::Id traND = static_cast<Dimensions::Transform::Id>(coord.ss().dimension()-1);

      auto pInVar = coord.ss().bwdPtr(traND);

      // Get the input data from hold
      if(edge.recoverInput())
      {
         coord.communicator().storage(traND).recoverBwd(pInVar);

      // Get the transfered input data
      } else
      {
         coord.communicator().receiveBackward(traND, pInVar);
      }

      // Get output storage
      auto pOutVar = coord.ss().fwdPtr(traND);
      auto pRecOutVar = coord.ss().fwdPtr(traND);
      if(edge.arithId() == Arithmetics::Set::id() || edge.arithId() == Arithmetics::SetNeg::id())
      {
         coord.communicator().storage(traND).recoverFwd(pOutVar);
      } else
      {
         coord.communicator().storage(traND).provideFwd(pOutVar);
         coord.communicator().storage(traND).recoverFwd(pRecOutVar);
      }

      ProfilerMacro_stop(Debug::Profiler::BWDNDIN);
      ProfilerMacro_start(Debug::Profiler::BWDNDTRA);

      // Compute projection transform for third dimension
      std::visit([&](auto&& pOut, auto&& pIn){coord.transformND().backward(pOut->rData(), pIn->data(), edge.opId());}, pOutVar, pInVar);

      ProfilerMacro_stop(Debug::Profiler::BWDNDTRA);
      ProfilerMacro_start(Debug::Profiler::BWDNDOUT);
      ProfilerMacro_start(Debug::Profiler::BWDNDOUTCOMM);

      // Hold temporary storage
      if(edge.holdInput())
      {
         coord.communicator().storage(traND).holdBwd(pInVar);

      // Free temporary input storage
      } else
      {
         coord.communicator().storage(traND).freeBwd(pInVar);
      }

      ProfilerMacro_stop(Debug::Profiler::BWDNDOUTCOMM);

      // Combine recovered output with new calculation
      if(std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
      {
         ProfilerMacro_start(Debug::Profiler::BWDNDOUTWORK);

         std::visit([&](auto&& pRecOut, auto&& pOut){Datatypes::FieldTools::combine(*pRecOut, *pOut, edge.arithId());}, pRecOutVar, pOutVar);

         ProfilerMacro_stop(Debug::Profiler::BWDNDOUTWORK);
         ProfilerMacro_start(Debug::Profiler::BWDNDOUTCOMM);

         coord.communicator().transferForward(traND, pRecOutVar);
         coord.communicator().storage(traND).freeFwd(pOutVar);

         ProfilerMacro_stop(Debug::Profiler::BWDNDOUTCOMM);
      } else
      {
         ProfilerMacro_start(Debug::Profiler::BWDNDOUTWORK);

         if(edge.arithId() == Arithmetics::SetNeg::id())
         {
            std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::BWDNDOUTWORK);
         ProfilerMacro_start(Debug::Profiler::BWDNDOUTCOMM);

         coord.communicator().transferForward(traND, pOutVar);

         ProfilerMacro_stop(Debug::Profiler::BWDNDOUTCOMM);
      }

      ProfilerMacro_stop(Debug::Profiler::BWDNDOUT);
      ProfilerMacro_stop(Debug::Profiler::BWDND);
   }

   void BackwardConfigurator::project1ND(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("BwdProject1ND");

      // Debugger message
      DebuggerMacro_msg("Project 1D", 4);

      ProfilerMacro_start(Debug::Profiler::BWD1D);
      ProfilerMacro_start(Debug::Profiler::BWD1DIN);

      // Get the input data from hold
      auto pInVar = coord.ss().bwdPtr(Dimensions::Transform::TRA1D);
      coord.communicator().storage(Dimensions::Transform::TRA1D).recoverBwd(pInVar);

      // Get output storage
      auto pOutVar = coord.ss().fwdPtr(Dimensions::Transform::TRA1D);
      auto pRecOutVar = coord.ss().fwdPtr(Dimensions::Transform::TRA1D);
      if(edge.recoverOutId() >= 0)
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).provideFwd(pOutVar);
         coord.communicator().storage(Dimensions::Transform::TRA1D).recoverFwd(pRecOutVar, edge.recoverOutId());
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).provideFwd(pOutVar);
      }

      ProfilerMacro_stop(Debug::Profiler::BWD1DIN);
      ProfilerMacro_start(Debug::Profiler::BWD1DTRA);

      // Compute projection transform for first dimension
      std::visit([&](auto&& pOut, auto&& pIn){coord.transform1D().backward(pOut->rData(), pIn->data(), edge.opId());}, pOutVar, pInVar);

      ProfilerMacro_stop(Debug::Profiler::BWD1DTRA);
      ProfilerMacro_start(Debug::Profiler::BWD1DOUT);
      ProfilerMacro_start(Debug::Profiler::BWD1DOUTCOMM);

      // Hold spectral input
      if(edge.holdInput())
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).holdBwd(pInVar);

      // Free spectral input
      } else
      {
         coord.communicator().storage(Dimensions::Transform::TRA1D).freeBwd(pInVar);
      }

      ProfilerMacro_stop(Debug::Profiler::BWD1DOUTCOMM);

      // Combine recovered output with new calculation
      if(std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
      {
         ProfilerMacro_start(Debug::Profiler::BWD1DOUTWORK);

         std::visit([&](auto&& pRecOut, auto&& pOut){Datatypes::FieldTools::combine(*pRecOut, *pOut, edge.arithId());}, pRecOutVar, pOutVar);

         ProfilerMacro_stop(Debug::Profiler::BWD1DOUTWORK);
         ProfilerMacro_start(Debug::Profiler::BWD1DOUTCOMM);

         coord.communicator().transferForward(Dimensions::Transform::TRA1D, pRecOutVar);
         coord.communicator().storage(Dimensions::Transform::TRA1D).freeFwd(pOutVar);

         ProfilerMacro_stop(Debug::Profiler::BWD1DOUTCOMM);
      } else
      {
         ProfilerMacro_start(Debug::Profiler::BWD1DOUTWORK);

         if(edge.arithId() == Arithmetics::SetNeg::id())
         {
            std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
         }

         ProfilerMacro_stop(Debug::Profiler::BWD1DOUTWORK);
         ProfilerMacro_start(Debug::Profiler::BWD1DOUTCOMM);

         coord.communicator().transferForward(Dimensions::Transform::TRA1D, pOutVar);

         ProfilerMacro_stop(Debug::Profiler::BWD1DOUTCOMM);
      }

      ProfilerMacro_stop(Debug::Profiler::BWD1DOUT);
      ProfilerMacro_stop(Debug::Profiler::BWD1D);
   }

}
}
