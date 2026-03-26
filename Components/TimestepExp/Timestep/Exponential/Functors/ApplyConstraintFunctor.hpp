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
      ApplyConstraintFunctor(const std::size_t t): timing(t) {};
      ~ApplyConstraintFunctor() = default;
      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt);
   private:
      std::size_t timing;
};

template <typename TEqIt>
void ApplyConstraintFunctor::operator()(const SpectralFieldId& myId, TEqIt& eqIt)
{
   // Apply constraint on solution
   eqIt->applyConstraint(myId.second, timing);
}


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
