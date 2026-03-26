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

      GetLinearInputFunctor(std::shared_ptr<TTsFunc> tsFunc, const std::size_t opId, const std::size_t regId, const std::size_t col, const Timestep::Interface::ScalarVariable_map& scalVar, const Timestep::Interface::VectorVariable_map& vectVar) : opId(opId), regId(regId), col(col), scalVar(scalVar), vectVar(vectVar), tsFunc(tsFunc){};
      ~GetLinearInputFunctor() = default;
      template <typename TEqIt>
      void operator()(ViewType tmpView, const SpectralFieldId& id, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i);
   private:
      const std::size_t opId;
      const std::size_t regId;
      const std::size_t col;
      const Timestep::Interface::ScalarVariable_map& scalVar;
      const Timestep::Interface::VectorVariable_map& vectVar;
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
template <typename TEqIt>
void GetLinearInputFunctor<TTsFunc>::operator()(ViewType tmpView, const SpectralFieldId& myId, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i)
{
   // Copy field values into timestepper input
   DecoupledZMatrix tmp(cinfo.tauN(i), cinfo.rhsCols(i));
   tmp.setZero();

   // Build range of operator
   auto r = make_range(cinfo.explicitRange(opId));

   // Loop over explicit fields
   for(auto& fIt: r)
   {
      DebuggerMacro_msg("Add " + ModelOperator::Coordinator::tag(opId) + " term from " + PhysicalNames::Coordinator::tag(fIt.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(fIt.second)) + ")", 7);

      // Get explicit input
      if(fIt.second == FieldComponents::Spectral::SCALAR)
      {
         std::visit(
               [&](auto&& p)
               {
               Equations::addExplicitTerm(*eqIt, opId, myId.second, tmp, 0, fIt, p->dom(0).perturbation(), i);
               }, scalVar.find(fIt.first)->second);
      } else
      {
         std::visit(
               [&](auto&& p)
               {
               Equations::addExplicitTerm(*eqIt, opId, myId.second, tmp, 0, fIt, p->dom(0).perturbation().comp(fIt.second), i);
               }, vectVar.find(fIt.first)->second);
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
