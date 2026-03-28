/**
 * @file TransferCorrectionFunctor.hpp
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
#include "Timestep/Exponential/CreateInfo.hpp"
#include "Timestep/Exponential/Functors/BaseFunctor.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Equations/CorrectSolution.hpp"
#include "QuICC/SolveTiming/After.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Transfer correction back to solver
 */
template <typename TTsFunc>
class TransferCorrectionFunctor
{
   public:
      TransferCorrectionFunctor(std::shared_ptr<TTsFunc> tsFunc, const std::size_t regId, const std::size_t col): regId(regId), col(col), tsFunc(tsFunc){};
      ~TransferCorrectionFunctor() = default;
      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const BaseFunctor::IdMap& idMap);
   protected:
      const std::size_t regId;
      const std::size_t col;
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
template <typename TEqIt>
void TransferCorrectionFunctor<TTsFunc>::operator()(const SpectralFieldId& myId, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const BaseFunctor::IdMap& idMap)
{
   // Apply constraint on solution
   auto changedSolution = eqIt->applyConstraint(myId.second, SolveTiming::After::id());

   // Update timestepper solver solution if constraint modified it
   if(changedSolution)
   {
      auto corr_ = eqIt->correctionConstraint(myId.second, SolveTiming::After::id());
      std::size_t matStart = 0;
      for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
      {
         auto info = createInfo(cinfo, i, idMap.at(myId));
         info.matStart = matStart;
         matStart += info.blockN;

         // Get effective corrections
         auto corr = Equations::correctSolution(eqIt->res(), cinfo, myId.second, corr_, i, 0);

         auto tsData = (*tsFunc)(info);
         auto&& pStepper = tsData.first;
         auto&& start = tsData.second;
         pStepper->correctData(corr, cinfo.galerkinN(i), cinfo.rhsCols(i), start, this->regId, this->col);
         pStepper->updateSolutions();
      }
   }
}


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
