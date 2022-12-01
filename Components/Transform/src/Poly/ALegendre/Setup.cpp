/**
 * @file Setup.cpp
 * @brief Source of polynomial transform setup class
 */

// System includes
//
#include "Kokkos.hpp"

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Setup.hpp"


namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

   Setup::Setup(const int size, const int specSize, const GridPurpose::Id purpose)
      : ::QuICC::Transform::Poly::Setup(size, specSize, purpose)
   {
      // Initialize fixtures
      ExternalLibrary::Kokkos::getInstance();
   }

   Setup::~Setup()
   {
   }

   void Setup::addIndex(const int slowIdx, const int mult)
   {
      ArrayI fastIdx = ArrayI::LinSpaced(this->specSize() - slowIdx, slowIdx, this->specSize() - 1);
      this->addIndex(slowIdx, mult, fastIdx);
   }

}
}
}
}
