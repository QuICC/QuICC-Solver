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
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/TransformConfigurators/BackwardConfigurator.hpp"

// Project includes
//
#include "QuICC/Arithmetics/SetNeg.hpp"
#include "QuICC/Arithmetics/None.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#ifdef QUICC_DEBUG
#include "QuICC/Transform/Backward/Coordinator.hpp"
#endif // QUICC_DEBUG
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   void BackwardConfigurator::prepareSpectral(const TransformTree& tree, Framework::Selector::VariantSharedScalarVariable& rScalar, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Transform::BackwardConfigurator::prepareSpectral");

      // Debugger message
      DebuggerMacro_msg("prepareSpectral (scalar)", 4);

      // Safety assert
      assert(tree.comp<FieldComponents::Spectral::Id>() == FieldComponents::Spectral::SCALAR);

      // Put scalar into temporary hold storage
      const auto traId = Dimensions::Transform::SPECTRAL;
      std::visit(
            [&](auto&& p)
            {
               coord.communicator().transferForward(traId, p->rDom(0).rTotal(), false);
            }, rScalar);
   }

   void BackwardConfigurator::prepareSpectral(const TransformTree& tree, Framework::Selector::VariantSharedVectorVariable& rVector, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("prepareSpectral (vector)", 4);

      Profiler::RegionFixture<2> fix("Transform::BackwardConfigurator::prepareSpectral");

      // Safety assert
      assert(tree.comp<FieldComponents::Spectral::Id>() != FieldComponents::Spectral::SCALAR);

      // Transfer scalar
      const auto traId = Dimensions::Transform::SPECTRAL;
      std::visit(
            [&](auto&& p)
            {
               coord.communicator().transferForward(traId, p->rDom(0).rTotal().rComp(tree.comp<FieldComponents::Spectral::Id>()), false);
            }, rVector);
   }

   void BackwardConfigurator::preparePhysical(const TransformTree&, const TransformTreeEdge& edge, Framework::Selector::VariantSharedScalarVariable& rScalar, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Transform::BackwardConfigurator::preparePhysical");

      // Debugger message
      DebuggerMacro_msg("preparePhysical (scalar)", 4);

      // Put scalar into temporary hold storage
      if(edge.fieldId() == FieldType::SCALAR)
      {
         std::visit([&](auto&& p){coord.communicator().holdPhysical(p->rDom(0).rPhys());}, rScalar);
      }
      // Put gradient component into temporary hold storage
      else if(edge.fieldId() == FieldType::GRADIENT)
      {
         std::visit([&](auto&& p){coord.communicator().holdPhysical(p->rDom(0).rGrad().rComp(edge.outId<FieldComponents::Physical::Id>()));}, rScalar);
      }
      // Put 2nd order gradient component into temporary hold storage
      else if(edge.fieldId() == FieldType::GRADIENT2)
      {
         std::visit([&](auto&& p){coord.communicator().holdPhysical(p->rDom(0).rGrad2().rComp(edge.outId<FieldComponents::Physical::Id>(0),edge.outId<FieldComponents::Physical::Id>(1)));}, rScalar);
      }
   }

   void BackwardConfigurator::preparePhysical(const TransformTree& tree, const TransformTreeEdge& edge, Framework::Selector::VariantSharedVectorVariable& rVector, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Transform::BackwardConfigurator::preparePhysical");

      // Debugger message
      DebuggerMacro_msg("preparePhysical (vector)", 4);

      // Put vector component into temporary hold storage
      if(edge.fieldId() == FieldType::VECTOR)
      {
         std::visit([&](auto&& p){coord.communicator().holdPhysical(p->rDom(0).rPhys().rComp(edge.outId<FieldComponents::Physical::Id>()));}, rVector);
      }
      // Put vector gradient component into temporary hold storage
      else if(edge.fieldId() == FieldType::GRADIENT)
      {
         std::visit([&](auto&& p){coord.communicator().holdPhysical(p->rDom(0).rGrad(tree.comp<FieldComponents::Spectral::Id>()).rComp(edge.outId<FieldComponents::Physical::Id>()));}, rVector);
      }
      // Put curl component into temporary hold storage
      else if(edge.fieldId() == FieldType::CURL)
      {
         std::visit([&](auto&& p){coord.communicator().holdPhysical(p->rDom(0).rCurl().rComp(edge.outId<FieldComponents::Physical::Id>()));}, rVector);
      }
   }

   void BackwardConfigurator::project1D(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Project 1D with operator (" + Backward::Coordinator::tag(edge.opId()) + ")", 6);

      const std::string profRegion = "Transform::BackwardConfigurator::project1D";
      const auto traId = Dimensions::Transform::TRA1D;
      const bool processOutput = true;
      BackwardConfigurator::genericProjection(edge, coord, traId, processOutput, profRegion);
   }

   void BackwardConfigurator::project2D(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Project 2D with operator (" + Backward::Coordinator::tag(edge.opId()) + ")", 6);

      const std::string profRegion = "Transform::BackwardConfigurator::project2D";
      const auto traId = Dimensions::Transform::TRA2D;
      const bool processOutput = true;
      BackwardConfigurator::genericProjection(edge, coord, traId, processOutput, profRegion);
   }

   void BackwardConfigurator::projectND(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Project ND with operator (" + Backward::Coordinator::tag(edge.opId()) + ")", 6);

      const std::string profRegion = "Transform::BackwardConfigurator::projectND";
      const auto traId = static_cast<Dimensions::Transform::Id>(coord.ss().dimension()-1);
      const bool processOutput = false;
      BackwardConfigurator::genericProjection(edge, coord, traId, processOutput, profRegion);
   }

   void BackwardConfigurator::genericProjection(const TransformTreeEdge& edge, TransformCoordinatorType& coord, const Dimensions::Transform::Id traId, const bool processOutput, const std::string profRegion)
   {
      Profiler::RegionFixture<2> fix(profRegion);
      Profiler::RegionStart<3> (profRegion + "-pre");

      auto pInVar = coord.ss().bwdPtr(traId);

      // Get the input data from hold
      if(edge.recoverInput())
      {
         coord.communicator().storage(traId).recoverBwd(pInVar);
      }
      // Get the transfered input data
      else
      {
         coord.communicator().receiveBackward(traId, pInVar);
      }

      // Get output storage
      auto pOutVar = coord.ss().fwdPtr(traId);
      auto pRecOutVar = coord.ss().fwdPtr(traId);
      if(processOutput)
      {
         if(edge.recoverOutId() >= 0)
         {
            coord.communicator().storage(traId).provideFwd(pOutVar);
            coord.communicator().storage(traId).recoverFwd(pRecOutVar, edge.recoverOutId());
         }
         else
         {
            coord.communicator().storage(traId).provideFwd(pOutVar);
         }
      }
      else
      {
         if(edge.arithId() == Arithmetics::Set::id() || edge.arithId() == Arithmetics::SetNeg::id())
         {
            coord.communicator().storage(traId).recoverFwd(pOutVar);
         }
         else
         {
            coord.communicator().storage(traId).provideFwd(pOutVar);
            coord.communicator().storage(traId).recoverFwd(pRecOutVar);
         }
      }

      Profiler::RegionStop<3> (profRegion + "-pre");
      Profiler::RegionStart<3> (profRegion + "-transform");

      // Compute projection transform
      std::visit([&](auto&& pOut, auto&& pIn){coord.transform(traId).backward(pOut->rData(), pIn->data(), edge.opId());}, pOutVar, pInVar);

      Profiler::RegionStop<3> (profRegion + "-transform");
      Profiler::RegionStart<3> (profRegion + "-post");
      Profiler::RegionStart<4> (profRegion + "-post" + "-comm_a");

      // Hold temporary storage
      if(edge.holdInput())
      {
         coord.communicator().storage(traId).holdBwd(pInVar);
      }
      // Free temporary storage
      else
      {
         coord.communicator().storage(traId).freeBwd(pInVar);
      }

      Profiler::RegionStop<4> (profRegion + "-post" + "-comm_a");

      if(processOutput)
      {
         // Combine recovered output with new calculation
         if(std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
         {
            Profiler::RegionStart<4> (profRegion + "-post" + "-work_b");

            std::visit([&](auto&& pRecOut, auto&& pOut){Datatypes::FieldTools::combine(*pRecOut, *pOut, edge.combinedArithId());}, pRecOutVar, pOutVar);

            Profiler::RegionStop<4> (profRegion + "-post" + "-work_b");
            Profiler::RegionStart<4> (profRegion + "-post" + "-comm_b");

            if(edge.combinedOutId() >= 0)
            {
               coord.communicator().storage(traId).holdFwd(pRecOutVar, edge.combinedOutId());
            }
            else
            {
               coord.communicator().transferForward(traId, pRecOutVar);
            }

            Profiler::RegionStop<4> (profRegion + "-post" + "-comm_b");
         }
         // Hold data for combination
         else if(edge.combinedOutId() >= 0)
         {
            Profiler::RegionStart<4> (profRegion + "-post" + "-work_c");

            if(edge.combinedArithId() == Arithmetics::SetNeg::id())
            {
               std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
            }

            Profiler::RegionStop<4> (profRegion + "-post" + "-work_c");
            Profiler::RegionStart<4> (profRegion + "-post" + "-comm_c");

            coord.communicator().storage(traId).holdFwd(pOutVar, edge.combinedOutId());

            Profiler::RegionStop<4> (profRegion + "-post" + "-comm_c");
         }

         // Transfer calculation
         if(edge.arithId() != Arithmetics::None::id())
         {
            Profiler::RegionStart<4> (profRegion + "-post" + "-work_d");

            if(edge.arithId() == Arithmetics::SetNeg::id())
            {
               std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
            }

            Profiler::RegionStop<4> (profRegion + "-post" + "-work_d");
            Profiler::RegionStart<4> (profRegion + "-post" + "-comm_d");

            coord.communicator().transferForward(traId, pOutVar);

            Profiler::RegionStop<4> (profRegion + "-post" + "-comm_d");
         }
         else if(edge.combinedOutId() < 0 || std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
         {
            Profiler::RegionStart<4> (profRegion + "-post" + "-comm_e");

            coord.communicator().storage(traId).freeFwd(pOutVar);

            Profiler::RegionStop<4> (profRegion + "-post" + "-comm_e");
         }
      }
      else
      {
         // Combine recovered output with new calculation
         if(std::visit([](auto&& pRecOut){return (pRecOut != 0);}, pRecOutVar))
         {
            Profiler::RegionStart<4> (profRegion + "-post" + "-work_f");

            std::visit([&](auto&& pRecOut, auto&& pOut){Datatypes::FieldTools::combine(*pRecOut, *pOut, edge.arithId());}, pRecOutVar, pOutVar);

            Profiler::RegionStop<4> (profRegion + "-post" + "-work_f");
            Profiler::RegionStart<4> (profRegion + "-post" + "-comm_f");

            coord.communicator().transferForward(traId, pRecOutVar);
            coord.communicator().storage(traId).freeFwd(pOutVar);

            Profiler::RegionStop<4> (profRegion + "-post" + "-comm_f");
         }
         else
         {
            Profiler::RegionStart<4> (profRegion + "-post" + "-work_g");

            if(edge.arithId() == Arithmetics::SetNeg::id())
            {
               std::visit([](auto&& pOut){Datatypes::FieldTools::negative(*pOut);}, pOutVar);
            }

            Profiler::RegionStop<4> (profRegion + "-post" + "-work_g");
            Profiler::RegionStart<4> (profRegion + "-post" + "-comm_g");

            coord.communicator().transferForward(traId, pOutVar);

            Profiler::RegionStop<4> (profRegion + "-post" + "-comm_g");
         }
      }

      Profiler::RegionStop<3> (profRegion + "-post");
   }

}
}
