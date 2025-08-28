/**
 * @file IReductor.cpp
 * @brief Source of the interface to a Finite Differences based projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/IReductor.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace FiniteDiff {

namespace Sphere {

namespace Reductor {

   IReductor::IReductor()
      : IOperator()
   {
      this->mProfileTag += "-Reductor";
   }

   MHDFloat IReductor::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += IOperator::requiredStorage();

      // Storage for the operators
      for(auto it = this->mOps.cbegin(); it != this->mOps.cend(); ++it)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(it->size());
      }

      // Storage for grid and weights
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(this->mGrid.size());
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}
}
}
}
