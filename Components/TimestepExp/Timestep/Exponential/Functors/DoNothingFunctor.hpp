/**
 * @file DoNothingFunctor.hpp
 * @brief
 */

#pragma once

// System includes
//

// Project includes
//

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Do Nothing
 */
class DoNothingFunctor
{
   public:
      DoNothingFunctor() = default;
      ~DoNothingFunctor() = default;
      template <typename ...Args>
      void operator()(const Args& ...){};
};


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
