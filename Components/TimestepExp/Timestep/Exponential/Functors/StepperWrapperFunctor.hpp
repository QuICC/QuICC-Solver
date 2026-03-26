/**
 * @file StepperWrapperFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <vector>

// Project includes
//
#include "Timestep/Exponential/TimestepperInfo.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Wrap current stepper
 */
template <typename TStepper>
class StepperWrapperFunctor
{
   public:
      typedef std::pair<TStepper*, std::size_t> ReturnType;

      StepperWrapperFunctor(const std::vector<std::size_t>& startArr): startArr(startArr), pStepper(nullptr){};
      ~StepperWrapperFunctor() = default;
      void setStepper(TStepper& ts);
      ReturnType operator()(const TimestepperInfo& info);
   private:
      std::vector<std::size_t> startArr;
      TStepper* pStepper;
};

template <typename TStepper>
typename StepperWrapperFunctor<TStepper>::ReturnType StepperWrapperFunctor<TStepper>::operator()(const TimestepperInfo& info)
{
   std::size_t start = startArr.at(info.fieldIndex) + info.matStart;
   auto ts = std::make_pair(pStepper, start);
   return ts;
}

template <typename TStepper>
void StepperWrapperFunctor<TStepper>::setStepper(TStepper& ts)
{
   this->pStepper = &ts;
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
