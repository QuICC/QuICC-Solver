/**
 * @file ProcessRangeFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/IteratorRange.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Process range of equations
 */
template <typename TBefore, typename TPrognostic, typename TAfter>
class ProcessRangeFunctor
{
   public:
      ProcessRangeFunctor(std::shared_ptr<TBefore> bFunc, std::shared_ptr<TPrognostic> pFunc, std::shared_ptr<TAfter> aFunc, const int fixedIt): _bFunc(bFunc), _pFunc(pFunc), _aFunc(aFunc), mRestricted(fixedIt != -1), mFixedIt(fixedIt) {};
      ~ProcessRangeFunctor() = default;
      template <typename TRange>
      void operator()(const TRange& eq_range);
   private:
      std::shared_ptr<TBefore> _bFunc;
      std::shared_ptr<TPrognostic> _pFunc;
      std::shared_ptr<TAfter> _aFunc;
      const bool mRestricted;
      const int mFixedIt;
};

template <typename TBefore, typename TPrognostic, typename TAfter>
template <typename TRange>
void ProcessRangeFunctor<TBefore,TPrognostic,TAfter>::operator()(const TRange& eq_range)
{
   TBefore& bFunc = *_bFunc;
   TPrognostic& pFunc = *_pFunc;
   TAfter& aFunc = *_aFunc;

   // Storage for information and identity
   SpectralFieldId myId;

   // Loop over equation range
   for(auto& eqIt: make_range(eq_range))
   {
      if(this->mRestricted && eqIt->options().it() != this->mFixedIt)
      {
         continue;
      }

      // Loop over spectral components
      for(auto& compId: make_range(eqIt->spectralRange()))
      {
         // Get field identity
         myId = std::make_pair(eqIt->name(), compId);

         // Process before prognostic equation
         bFunc(myId, eqIt);

         // Process prognostic equation
         if(eqIt->solveTiming() == SolveTiming::Prognostic::id())
         {
            pFunc(myId, eqIt);
         }

         // Process after prognostic equation
         aFunc(myId, eqIt);
      }
   }
}


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
