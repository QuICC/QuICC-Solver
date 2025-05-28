/**
 * @file Setup.cpp
 * @brief Source of finite differences transform setup class
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/FiniteDiff/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace FiniteDiff {

   Setup::Setup(const int size, const GridPurpose::Id purpose)
      : TransformSetup(size, size, purpose)
   {
   }

} // namespace FiniteDiff
} // namespace Transform
} // namespace QuICC
