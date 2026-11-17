/**
 * @file IFftwBackend.cpp
 * @brief Source of the interface for a generic FFTW based backend
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Backend/Fftw/IFftwBackend.hpp"
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

   void IFftwBackend::applyFft(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in) const
   {
      fftw_execute_r2r(this->mPlan, const_cast<MHDFloat *>(in.data()), rOut.data());
   }

   void IFftwBackend::applyFft(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in) const
   {
      fftw_execute_dft_c2r(this->mPlan, reinterpret_cast<fftw_complex* >(const_cast<MHDComplex *>(in.data())), rOut.data());
   };

   void IFftwBackend::applyFft(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& in) const
   {
      fftw_execute_dft_r2c(this->mPlan, const_cast<MHDFloat*>(in.data()), reinterpret_cast<fftw_complex* >(rOut.data()));
   };

   void IFftwBackend::applyFft(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in) const
   {
      fftw_execute_dft(this->mPlan, reinterpret_cast<fftw_complex *>(const_cast<MHDComplex*>(in.data())), reinterpret_cast<fftw_complex *>(rOut.data()));
   };

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

   void IFftwBackend::applyFft() const {std::logic_error("Backend not implemented.");};

}
}
}
}
}
