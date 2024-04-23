/**
 * @file IBesselReductor.cpp
 * @brief Source of the interface to a spherical Bessel based reduction operator (e.g. energy)
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Reductor/IBesselReductor.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   IBesselReductor::IBesselReductor()
      : IBesselOperator()
   {
      this->mProfileTag += "-Reductor";
   }

   MHDFloat IBesselReductor::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += IBesselOperator::requiredStorage();

      // Storage for the operators
      for(auto it = this->mOps.cbegin(); it != this->mOps.cend(); ++it)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(it->size());
      }

      // Storage for grid and weights
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(this->mGrid.size());
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(this->mWeights.size());
#endif // QUICC_STORAGEPROFILE

      return mem;
   }
}
}
}
}
}
