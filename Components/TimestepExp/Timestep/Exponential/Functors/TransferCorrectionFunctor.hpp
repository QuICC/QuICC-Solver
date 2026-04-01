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
#include "Timestep/Exponential/Functors/FunctorData.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Equations/CorrectSolution.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#endif //QUICC_DEBUG

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
      TransferCorrectionFunctor(std::shared_ptr<FunctorData> spData, std::shared_ptr<TTsFunc> tsFunc, const std::size_t regId, const std::size_t col): regId(regId), col(col), spData(spData), tsFunc(tsFunc){};
      ~TransferCorrectionFunctor() = default;
      void operator()(const SpectralFieldId& id, const Equations::CouplingInformation& cinfo, const BaseFunctor::IdMap& idMap);
   protected:
      const std::size_t regId;
      const std::size_t col;
      std::shared_ptr<FunctorData> spData;
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
void TransferCorrectionFunctor<TTsFunc>::operator()(const SpectralFieldId& myId, const Equations::CouplingInformation& cinfo, const BaseFunctor::IdMap& idMap)
{
   const auto& data = *spData;

   const auto timeId = SolveTiming::After::id();
   bool changedSolution = false;

   // Apply constraint on solution
   if(data.constraints.count(myId) > 0)
   {
      DebuggerMacro_msg("Apply constraint kernel for " + PhysicalNames::Coordinator::tag(myId.first) + "(" + QuICC::Tools::IdToHuman::toString(myId.second) + ") at " + SolveTiming::Coordinator::tag(timeId) , 6);

      changedSolution = true;
      data.constraints.at(myId)->apply(timeId);
   }

   // Update timestepper solver solution if constraint modified it
   if(changedSolution)
   {
      std::vector<std::tuple<MHDVariant,int,int,int>> corr_;

      if(data.constraints.count(myId) > 0)
      {
         DebuggerMacro_msg("Correction from constraint kernel for " + PhysicalNames::Coordinator::tag(myId.first) + "(" + QuICC::Tools::IdToHuman::toString(myId.second) + ") at " + SolveTiming::Coordinator::tag(timeId) , 6);

         corr_ = data.constraints.at(myId)->correction(timeId);
      }

      std::size_t matStart = 0;
      for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
      {
         auto info = createInfo(cinfo, i, idMap.at(myId));
         info.matStart = matStart;
         matStart += info.blockN;

         // Get effective corrections
         auto corr = Equations::correctSolution(data.res(), cinfo, corr_, i, 0);

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
