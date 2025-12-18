/**
 * @file ILinearMapSpectrum.cpp
 * @brief Source of the interface for a generic FFT based Chebyshev spectrum reductor
 */

// System includes
//
#include <cassert>
#include <stdexcept>


// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/ILinearMapSpectrum.hpp"


namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   void ILinearMapSpectrum::initBackend() const
   {
      this->mBackend.init(*this->mspSetup);
   }

   void ILinearMapSpectrum::transform(Matrix& rOut, const MatrixZ& in) const
   {
      rOut = in.array().abs2();
   }

   // anelastic version:
   void ILinearMapSpectrum::transform(Matrix& rOut, const MatrixZ& in, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
   {      
      // Energy version
      //assert(this->isInitialized());
      //assert(rOut.cols() == this->outCols()); // fails
      //assert(rOut.rows() == this->outRows());

      auto eGrid = this->mBackend.getEGrid();
      
      auto rho = pF->evaluateLP(eGrid,0,0);

      auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
      auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
      auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);

      // real part:
      this->mBackend.input(tmpIn,in,true); 
      // field in physical space
      this->mBackend.applyFft(tmpOut, tmpIn);  
      // divide by sqrt(rho)
      tmpOut = tmpOut.array().colwise() / rho.array().pow(0.5);
      // calculate spectra
      this->mBackend.applyFwdFft(tmpOut, tmpOut); 
      // adjust for fft scaling:
      tmpOut = tmpOut.array()*(this->mBackend.getFftScaling());
      // square it
      this->mBackend.square(tmpSquare, tmpOut, true); // tmpSquare now contains the spectral coefficients, real part


      // imaginary part:
      this->mBackend.input(tmpIn,in,false);
      // field in physical space
      this->mBackend.applyFft(tmpOut, tmpIn); 
      // divide by sqrt(rho)
      tmpOut = tmpOut.array().colwise() / rho.array().pow(0.5);
      // calculate spectra
      this->mBackend.applyFwdFft(tmpOut, tmpOut); 
      // adjust for fft scaling:
      tmpOut = tmpOut.array()*(this->mBackend.getFftScaling());
      // square it
      this->mBackend.square(tmpSquare, tmpOut, false); // tmpSquare now contains the spectral coefficients, imaginary and real part

      rOut = tmpSquare;
   }


   void ILinearMapSpectrum::transform(Matrix& rOut, const Matrix& in) const
   {
      rOut = in.array().abs2();
   }

   void ILinearMapSpectrum::transform(Matrix& rOut, const Matrix& in, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
   {
      throw std::logic_error("Anelastic spectrum not yet implemented");
   }
   
   void ILinearMapSpectrum::transform(MatrixZ&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Chebyshev FFT energy reductor");
   }

   void ILinearMapSpectrum::transform(MatrixZ&, const Matrix&) const
   {
      throw std::logic_error("Data is not compatible with Chebyshev FFT energy reductor");
   }
   
   MHDFloat ILinearMapSpectrum::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   int ILinearMapSpectrum::outRows() const
   {
      return this->mspSetup->fwdSize();
   }

   int ILinearMapSpectrum::outCols() const
   {
      return this->mspSetup->blockSize();
   }

}
}
}
}
}
}
