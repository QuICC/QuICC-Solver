/**
 * @file ApplyConstraitFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "Timestep/Exponential/Functors/FunctorData.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/SolveTiming/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#endif //QUICC_DEBUG

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Apply constraint functor
 */
class ApplyConstraintFunctor
{
   public:
      ApplyConstraintFunctor(std::shared_ptr<FunctorData> spData, const std::size_t t): spData(spData), timeId(t) {};
      ~ApplyConstraintFunctor() = default;
      void operator()(const SpectralFieldId& id);
   private:
      std::shared_ptr<FunctorData> spData;
      std::size_t timeId;
};

inline void ApplyConstraintFunctor::operator()(const SpectralFieldId& myId)
{
   const auto& eqData = *spData;

   // Apply constraint on solution
   if(eqData.constraints.count(myId) > 0)
   {
      DebuggerMacro_msg("Apply constraint kernel for " + PhysicalNames::Coordinator::tag(myId.first) + "(" + QuICC::Tools::IdToHuman::toString(myId.second) + ") at " + SolveTiming::Coordinator::tag(timeId) , 6);

      eqData.constraints.at(myId)->apply(timeId);
   }
}


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
