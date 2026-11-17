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
#include "Timestep/Exponential/Functors/FunctorData.hpp"

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

      InitSolutionFunctor(std::shared_ptr<FunctorData> pData, std::shared_ptr<TTsFunc> tsFunc): pData(pData), tsFunc(tsFunc){};
      ~InitSolutionFunctor() = default;
      void operator()(ViewType tmpView, const SpectralFieldId& id, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i);
   protected:
      std::shared_ptr<FunctorData> pData;
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
void InitSolutionFunctor<TTsFunc>::operator()(ViewType tmpView, const SpectralFieldId& myId, const Equations::CouplingInformation& , const TimestepperInfo& info, const std::size_t i)
{
   const auto& data = *pData;
   const auto& cinfo = data.cInfos.at(myId);

   if(cinfo.isGalerkin())
   {
      assert(data.fields.count(myId) > 0);
      Equations::solveStencilUnknown(data.res(), data.cInfos.at(myId), myId, *data.fields.at(myId), tmpView, i, 0, data.backend(), data.bcIdMap, data.eqParamsMap);
   }
   else
   {
      assert(data.fields.count(myId) > 0);
      Equations::copyUnknown(data.res(), data.cInfos.at(myId), *data.fields.at(myId), tmpView, i, 0, true, true, true);
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
