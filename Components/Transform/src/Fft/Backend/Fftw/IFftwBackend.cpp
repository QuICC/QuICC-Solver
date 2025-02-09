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
      QuICC::Fft::Fftw::Library::getInstance();
   }

   void IFftwBackend::applyFft(Matrix& phys, const Matrix& mods) const
   {
      fftw_execute_r2r(this->mPlan, const_cast<MHDFloat *>(mods.data()), phys.data());
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

   // to be removed
   void IFftwBackend::applyFft() const {std::logic_error("Backend not implemented.");};
   void IFftwBackend::applyFft(Matrix&, const MatrixZ&) const {std::logic_error("Backend not implemented.");};
   void IFftwBackend::applyFft(MatrixZ&, const Matrix&) const {std::logic_error("Backend not implemented.");};
   void IFftwBackend::applyFft(MatrixZ&, const MatrixZ&) const {std::logic_error("Backend not implemented.");};

}
}
}
}
}
