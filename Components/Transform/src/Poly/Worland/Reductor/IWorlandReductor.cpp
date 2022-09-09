/**
 * @file IWorlandReductor.cpp
 * @brief Source of the interface to a Worland based reduction operator (e.g. energy)
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandReductor.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   IWorlandReductor::IWorlandReductor()
      : IWorlandOperator()
   {
      this->mProfileTag += "-Reductor";
   }

   IWorlandReductor::~IWorlandReductor()
   {
   }

   MHDFloat IWorlandReductor::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += IWorlandOperator::requiredStorage();

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
