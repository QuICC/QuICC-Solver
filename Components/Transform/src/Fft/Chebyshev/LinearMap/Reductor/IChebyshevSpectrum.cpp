/**
 * @file IChebyshevSpectrum.cpp
 * @brief Source of the interface for a generic FFT based Chebyshev spectrum reductor
 */

// System includes
//
#include <cassert>
#include <stdexcept>
// ****************
//Stuff that needs to be removed later
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <iostream>
// ****************

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/IChebyshevSpectrum.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   IChebyshevSpectrum::IChebyshevSpectrum()
   {
   }

   IChebyshevSpectrum::~IChebyshevSpectrum()
   {
   }

   void IChebyshevSpectrum::initBackend() const
   {
      this->mBackend.init(*this->mspSetup);
   }

   void IChebyshevSpectrum::transform(Matrix& rOut, const MatrixZ& in) const
   {
      rOut = in.array().abs2();

      //std::cerr << "in = \n";
      //std::cerr << in << "\n";
      //std::cerr << "\n";
      //std::cerr << "rOut = \n";
      //std::cerr << rOut << "\n";
      //std::cerr << "\n";

      /*
      // Energy version
      assert(this->isInitialized());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() == this->outRows());

      auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
      auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
      auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);
      this->applyPreOperator(tmpIn, in, true);
      this->mBackend.applyFft(tmpOut, tmpIn);
      this->mBackend.square(tmpSquare, tmpOut, true);
      this->applyPreOperator(tmpIn, in, false);
      this->mBackend.applyFft(tmpOut, tmpIn);
      this->mBackend.square(tmpSquare, tmpOut, false);
      this->mBackend.applyFwdFft(tmpOut, tmpSquare);
      this->applyPostOperator(rOut, tmpOut);
      */
   }

   // anelastic version:
   void IChebyshevSpectrum::transform(Matrix& rOut, const MatrixZ& in, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
   {
      //throw std::logic_error("Anelastic spectrum not yet implemented");
      
      // Energy version
      //assert(this->isInitialized());
      //assert(rOut.cols() == this->outCols()); // fails
      //assert(rOut.rows() == this->outRows());

      auto eGrid = this->mBackend.getEGrid();
      
      auto rho = pF->evaluate(eGrid,0,0);

      auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
      //auto& tmpOutRe = this->mBackend.getStorage(StorageKind::out);
      auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
      auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);

      // real part:
      //this->applyPreOperator(tmpIn, in, true);
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
      // store real part of spectra in a new temporary matrix:
      //Matrix tmpOutIm = tmpOut;

      // imaginary part:
      //this->applyPreOperator(tmpIn, in, false);
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
      //this->mBackend.square(tmpSquare, tmpOut, false);

      //std::cerr << "rOut = \n";
      //std::cerr << rOut << "\n";
      //std::cerr << "\n";

      //this->applyPostOperator(rOut, tmpOut);  // I don't think we need this for Chebyshev spectra calculation
      // HOWEVER::::
      // the result of FwdFftw is 2*Ngrid too big. Postoperator does something to fix this, I guess.
   }


   void IChebyshevSpectrum::transform(Matrix& rOut, const Matrix& in) const
   {
      rOut = in.array().abs2();
      /*
      assert(this->isInitialized());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() == this->outRows());

      auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
      auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
      auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);
      this->applyPreOperator(tmpIn, in);
      this->mBackend.applyFft(tmpOut, tmpIn);
      this->mBackend.square(tmpSquare, tmpOut, true);
      this->mBackend.applyFwdFft(tmpOut, tmpSquare);
      this->applyPostOperator(rOut, tmpOut);
      */
   }

   void IChebyshevSpectrum::transform(Matrix& rOut, const Matrix& in, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
   {
      throw std::logic_error("Anelastic spectrum not yet implemented");
      /*
      // Energy version
      assert(this->isInitialized());
      //assert(rOut.cols() == this->outCols());
      //assert(rOut.rows() == this->outRows());

      auto eGrid = this->mBackend.getEGrid();
      
      auto rho = pF->evaluate(eGrid,0,0);

      auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
      auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
      auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);
      this->applyPreOperator(tmpIn, in);
      this->mBackend.applyFft(tmpOut, tmpIn);
      this->mBackend.square(tmpSquare, tmpOut, true);
      //tmpSquare = tmpSquare.array().colwise() / rho.array(); // divides energy by rho
      this->mBackend.applyFwdFft(tmpOut, tmpSquare);
      this->applyPostOperator(rOut, tmpOut);
      */
   }

   void IChebyshevSpectrum::transform(MatrixZ&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Chebyshev FFT energy reductor");
   }

   void IChebyshevSpectrum::transform(MatrixZ&, const Matrix&) const
   {
      throw std::logic_error("Data is not compatible with Chebyshev FFT energy reductor");
   }

   MHDFloat IChebyshevSpectrum::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   int IChebyshevSpectrum::outRows() const
   {
      return this->mspSetup->blockSize();
   }

   int IChebyshevSpectrum::outCols() const
   {
      return 1;
   }

}
}
}
}
}
}
