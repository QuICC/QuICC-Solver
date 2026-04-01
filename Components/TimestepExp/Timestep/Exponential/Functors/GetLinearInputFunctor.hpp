/**
 * @file GetLinearInputFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <memory>

// Project includes
//
#include "Memory/MemoryResource.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "Profiler/Interface.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Equations/ExplicitTerm.hpp"
#include "QuICC/Timestep/Interface.hpp"
#include "Timestep/Exponential/TimestepperInfo.hpp"
#include "Timestep/Exponential/details/TimesteppperTools.hpp"
#include "Timestep/Exponential/Functors/FunctorData.hpp"
#include "View/Attributes.hpp"
#include "View/ViewDense.hpp"
#include "QuICC/IteratorRange.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Get linear input functor
 */
template <typename TTsFunc>
class GetLinearInputFunctor
{
   public:
      using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
      typedef View::View<MHDComplex, View::Attributes<dense2D>> ViewType;

      GetLinearInputFunctor(std::shared_ptr<FunctorData> pData, std::shared_ptr<TTsFunc> tsFunc, const std::size_t opId, const std::size_t regId, const std::size_t col, const Timestep::Interface::ScalarVariable_map& scalVar, const Timestep::Interface::VectorVariable_map& vectVar) : opId(opId), regId(regId), col(col), scalVar(scalVar), vectVar(vectVar), pData(pData), tsFunc(tsFunc){};
      ~GetLinearInputFunctor() = default;
      void operator()(ViewType tmpView, const SpectralFieldId& id, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i);
   private:
      const std::size_t opId;
      const std::size_t regId;
      const std::size_t col;
      const Timestep::Interface::ScalarVariable_map& scalVar;
      const Timestep::Interface::VectorVariable_map& vectVar;
      std::shared_ptr<FunctorData> pData;
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
void GetLinearInputFunctor<TTsFunc>::operator()(ViewType tmpView, const SpectralFieldId& myId, const Equations::CouplingInformation&, const TimestepperInfo& info, const std::size_t i)
{
   const auto& data = *pData;
   assert(data.cInfos.count(myId) > 0);
   const auto& cinfo = data.cInfos.at(myId);

   // Copy field values into timestepper input
   DecoupledZMatrix tmp(cinfo.tauN(i), cinfo.rhsCols(i));
   tmp.setZero();

   if(data.exDTerm.count(opId) > 0 && data.exDTerm.at(opId).count(myId) > 0)
   {
      for(auto&& [exId, mats]: data.exDTerm.at(opId).at(myId))
      {
         DebuggerMacro_msg("Add real " + ModelOperator::Coordinator::tag(opId) + " term from " + PhysicalNames::Coordinator::tag(exId.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(exId.second)) + ")", 7);

         Equations::addExplicitTerm(data.res(), cinfo, mats.at(i), tmp, 0, *data.fields.at(exId), i);
      }
   }

   if(data.exZTerm.count(opId) > 0 && data.exZTerm.at(opId).count(myId) > 0)
   {
      for(auto&& [exId, mats]: data.exZTerm.at(opId).at(myId))
      {
         DebuggerMacro_msg("Add complex " + ModelOperator::Coordinator::tag(opId) + " term from " + PhysicalNames::Coordinator::tag(exId.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(exId.second)) + ")", 7);

         Equations::addExplicitTerm(data.res(), cinfo, mats.at(i), tmp, 0, *data.fields.at(exId), i);
      }
   }

   details::computeSet(tmpView, tmp, cinfo.tauN(i) - cinfo.galerkinN(i));

   auto tsData = (*tsFunc)(info);
   auto&& pStepper = tsData.first;
   auto&& start = tsData.second;
   pStepper->addData(tmpView, start, this->regId, this->col);
}


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
