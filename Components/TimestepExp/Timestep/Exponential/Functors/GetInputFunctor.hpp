/**
 * @file GetInputFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <memory>

// Project includes
//
#include "Memory/MemoryResource.hpp"
#include "View/Attributes.hpp"
#include "View/ViewDense.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Equations/CopyUnknown.hpp"
#include "QuICC/Equations/AddSource.hpp"
#include "Timestep/Exponential/TimestepperInfo.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

template <typename TTsFunc>
class GetInputFunctor
{
   public:
      using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
      typedef View::View<MHDComplex, View::Attributes<dense2D>> ViewType;

      GetInputFunctor(std::shared_ptr<TTsFunc> tsFunc, const std::size_t regId, const std::size_t col): regId(regId), col(col), tsFunc(tsFunc){};
      ~GetInputFunctor() = default;
      template <typename TEqIt>
      void operator()(ViewType tmpView, const SpectralFieldId& id, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i);

   protected:
      const std::size_t regId;
      const std::size_t col;
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
template <typename TEqIt>
void GetInputFunctor<TTsFunc>::operator()(ViewType tmpView, const SpectralFieldId& myId, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i)
{
   // Copy field values into timestepper input
   if(cinfo.hasNonlinear())
   {
      std::visit(
            [&](auto&& p)
            {
            Equations::copyUnknown(eqIt->res(), eqIt->couplingInfo(myId.second), p->dom(0).perturbation(), myId.second, tmpView, i, 0, true, true, false);
            }, eqIt->spUnknown());
   }

   // Add source term
   if(cinfo.hasSource())
   {
      std::visit(
            [&](auto&& p)
            {
            Equations::addSource(eqIt->res(), cinfo, eqIt->sourceKernel(myId.second), p->dom(0).perturbation(), myId.second, tmpView, i, 0);
            }, eqIt->spUnknown());
   }

   // Add value to RHS
   auto tsData = (*tsFunc)(info);
   auto&& pStepper = tsData.first;
   auto&& start = tsData.second;
   pStepper->addQiData(tmpView, start, this->regId, this->col);

   // Enforce BC
   const auto& rows = tmpView.dims()[0];
   const auto& cols = tmpView.dims()[1];
   pStepper->enforceBoundaryConditions(rows, cols, start, this->regId, this->col);

   // If required set inhomogenous boundary condition value
   if(cinfo.hasBoundaryValue())
   {
      throw std::logic_error("Inhomogeneous boundary conditions are not supported!");
   }
}


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
