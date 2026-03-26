/**
 * @file InitSolutionFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Equations/SolveStencilUnknown.hpp"
#include "Memory/MemoryResource.hpp"
#include "View/Attributes.hpp"
#include "View/ViewDense.hpp"
#include "Timestep/Exponential/TimestepperInfo.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Init solution functor
 */
template <typename TTsFunc>
class InitSolutionFunctor
{
   public:
      using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
      typedef View::View<MHDComplex, View::Attributes<dense2D>> ViewType;

      InitSolutionFunctor(std::shared_ptr<TTsFunc> tsFunc): tsFunc(tsFunc){};
      ~InitSolutionFunctor() = default;
      template <typename TEqIt>
      void operator()(ViewType tmpView, const SpectralFieldId& id, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i);
   protected:
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
template <typename TEqIt>
void InitSolutionFunctor<TTsFunc>::operator()(ViewType tmpView, const SpectralFieldId& myId, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i)
{
   if(cinfo.isGalerkin())
   {
      Equations::solveStencilUnknown(*eqIt, myId.second, tmpView, i, 0);
   }
   else
   {
      std::visit(
            [&](auto&& p)
            {
            Equations::copyUnknown(*eqIt, p->dom(0).perturbation(), myId.second, tmpView, i, 0, true, true, true);
            }, eqIt->spUnknown());
   }

   auto tsData = (*tsFunc)(info);
   auto&& pStepper = tsData.first;
   auto&& start = tsData.second;
   pStepper->setSolution(tmpView, start);
   pStepper->updateSolutions();
}


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
