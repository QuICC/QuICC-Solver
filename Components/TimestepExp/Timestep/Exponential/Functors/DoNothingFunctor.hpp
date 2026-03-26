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
      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt){};
      
};


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
