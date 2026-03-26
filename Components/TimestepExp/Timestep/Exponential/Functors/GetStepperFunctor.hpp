/**
 * @file GetStepperFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//

// Project includes
//
#include "Timestep/Exponential/TimestepperInfo.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Get current stepper
 */
template <typename TCoord>
class GetStepperFunctor
{
   public:
      typedef std::pair<typename TCoord::TimestepperType*, std::size_t> ReturnType;

      GetStepperFunctor(TCoord& coord): coord(coord) {};
      ~GetStepperFunctor() = default;
      ReturnType operator()(const TimestepperInfo& info);
   private:
      TCoord& coord;
};

template <typename TCoord>
typename GetStepperFunctor<TCoord>::ReturnType GetStepperFunctor<TCoord>::operator()(const TimestepperInfo& info)
{
   return coord.getStepper(info);
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
