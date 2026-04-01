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
#include "QuICC/Equations/StoreSolution.hpp"
#include "Timestep/Exponential/TimestepperInfo.hpp"
#include "View/Attributes.hpp"
#include "View/ViewDense.hpp"
#include "Timestep/Exponential/details/TimesteppperTools.hpp"
#include "Timestep/Exponential/Functors/FunctorData.hpp"

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

      /**
       * @brief ctor
       */
      TransferOutputFunctor(std::shared_ptr<FunctorData> spData, std::shared_ptr<TTsFunc> tsFunc, const std::size_t regId, const std::size_t col) : regId(regId), col(col), spData(spData), tsFunc(tsFunc){};

      /**
       * @brief dtor
       */
      ~TransferOutputFunctor() = default;

      void operator()(ViewType tmpView, const SpectralFieldId& id, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i);

   protected:
      const std::size_t regId;
      const std::size_t col;
      std::shared_ptr<FunctorData> spData;
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
void TransferOutputFunctor<TTsFunc>::operator()(ViewType tmpView, const SpectralFieldId& myId, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i)
{
   auto& eqData = *spData;
   auto tsData = (*tsFunc)(info);
   auto&& pStepper = tsData.first;
   auto&& start = tsData.second;
   pStepper->getData(tmpView, start, this->regId, this->col);

   DecoupledZMatrix tmp(cinfo.galerkinN(i), cinfo.rhsCols(i));
   details::computeSet(tmp, tmpView, 0);

   const SparseMatrix* pOp;
   if(cinfo.isGalerkin())
   {
      assert(eqData.stencils.count(myId) > 0);
      pOp = &eqData.stencils.at(myId).at(i);
   }
   assert(eqData.fields.count(myId) > 0);
   assert(eqData.solups.count(myId) > 0);
   Equations::storeSolution(*eqData.fields.at(myId), eqData.res(), cinfo, pOp, eqData.solups.at(myId), tmp, i, 0);
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
