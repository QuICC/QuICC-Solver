/**
 * @file Setup.cpp
 * @brief Source of polynomial transform setup class
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Setup.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Poly {

   Setup::Setup(const int size, const int specSize, const GridPurpose::Id purpose)
      : TransformSetup(size, specSize, purpose)
   {
   }

   Setup::~Setup()
   {
   }

}
}
}
