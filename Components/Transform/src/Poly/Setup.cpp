/**
 * @file Setup.cpp
 * @brief Source of polynomial transform setup class
 */

// System includes
//
#include <stdexcept>
#include "Kokkos.hpp"

// Project includes
//
#include "QuICC/Transform/Poly/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

   Setup::Setup(const int size, const int specSize, const GridPurpose::Id purpose)
      : TransformSetup(size, specSize, purpose)
   {
       // Initialize fixtures
       ExternalLibrary::Kokkos::getInstance();
   }

}
}
}
