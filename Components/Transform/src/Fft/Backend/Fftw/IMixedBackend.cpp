/**
 * @file IMixedBackend.cpp
 * @brief Source of the interface for a generic FFTW based mixed backend
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/IMixedBackend.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   IMixedBackend::IMixedBackend()
   {
   }

   IMixedBackend::~IMixedBackend()
   {
   }

   void IMixedBackend::init(const SetupType& setup) const
   {
      this->mSpecSize = setup.specSize();
   }

   Array IMixedBackend::positiveK() const
   {
      return Array::LinSpaced(this->mSpecSize, 0, this->mSpecSize-1);
   }

   MHDFloat IMixedBackend::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += IFftwBackend::requiredStorage();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDComplex>::BYTES*this->mTmp.size());
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   MatrixZ& IMixedBackend::getStorage() const
   {
      return this->mTmp;
   }

}
}
}
}
}
