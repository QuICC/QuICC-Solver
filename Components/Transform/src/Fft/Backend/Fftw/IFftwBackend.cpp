/**
 * @file IFftwBackend.cpp
 * @brief Source of the interface for a generic FFTW based backend
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/IFftwBackend.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   IFftwBackend::IFftwBackend()
      : mPlan(NULL)
   {
      this->initLibrary();
   }

   IFftwBackend::~IFftwBackend()
   {
      this->cleanupFft();
   }

   void IFftwBackend::initLibrary() const
   {
      // FFTW Fixture
      Library::getInstance();
   }

   void IFftwBackend::cleanupFft()
   {
      // Destroy plan
      if(this->mPlan)
      {
         fftw_destroy_plan(this->mPlan);
      }

   }

   MHDFloat IFftwBackend::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}
}
}
}
