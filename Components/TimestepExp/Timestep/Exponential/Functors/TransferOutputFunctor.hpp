/**
 * @file TransferOutputFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <memory>

// Project includes
//
#include "Memory/MemoryResource.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "Timestep/Exponential/TimestepperInfo.hpp"
#include "View/Attributes.hpp"
#include "View/ViewDense.hpp"
#include "Timestep/Exponential/details/TimesteppperTools.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Transfer output functor
 */
template <typename TTsFunc>
class TransferOutputFunctor
{
   public:
      using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
      typedef View::View<MHDComplex, View::Attributes<dense2D>> ViewType;

      TransferOutputFunctor(std::shared_ptr<TTsFunc> tsFunc, const std::size_t regId, const std::size_t col) : regId(regId), col(col), tsFunc(tsFunc){};
      ~TransferOutputFunctor() = default;
      template <typename TEqIt>
      void operator()(ViewType tmpView, const SpectralFieldId& id, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i);

   protected:
      const std::size_t regId;
      const std::size_t col;
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
template <typename TEqIt>
void TransferOutputFunctor<TTsFunc>::operator()(ViewType tmpView, const SpectralFieldId& myId, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i)
{
   auto tsData = (*tsFunc)(info);
   auto&& pStepper = tsData.first;
   auto&& start = tsData.second;
   pStepper->getData(tmpView, start, this->regId, this->col);

   DecoupledZMatrix tmp(cinfo.galerkinN(i), cinfo.rhsCols(i));
   details::computeSet(tmp, tmpView, 0);

   eqIt->storeSolution(myId.second, tmp, i, 0);
}


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
