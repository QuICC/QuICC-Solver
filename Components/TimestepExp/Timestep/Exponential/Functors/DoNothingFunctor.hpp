/**
 * @file DoNothingFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"

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
      void operator()(const SpectralFieldId& id){};
      template <typename TEq> 
      void operator()(const SpectralFieldId& id, const TEq&){};
};


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
