/**
 * @file IChebyshevEnergy.cpp
 * @brief Source of the interface for a generic FFT based Chebyshev energy reductor
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//
// ****************
//Stuff that needs to be removed later
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <iostream>
// ****************

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/IChebyshevEnergy.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   IChebyshevEnergy::IChebyshevEnergy()
   {
   }

   IChebyshevEnergy::~IChebyshevEnergy()
   {
   }

   void IChebyshevEnergy::initBackend() const
   {
      this->mBackend.init(*this->mspSetup);
   }

   // anelastic version (overload not possible without having to alter a bunch of other classes)
   void IChebyshevEnergy::initBackendAnelastic(std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
   {
      // Set the extra size BEFORE backend initialization
      this->mExtraSize=pF->nN();
      // set it in the backend too
      this->mBackend.setExtraSize(pF->nN()); 

      this->mBackend.init(*this->mspSetup, pF);

   }

   int IChebyshevEnergy::getExtraSize() const
   {
      return this->mExtraSize;
   }

   void IChebyshevEnergy::transform(Matrix& rOut, const MatrixZ& in) const
   {
      // NB: fwd is pysical to spectral
      // input and output sizes = (setup.bwdSize rows) x (blockSize cols)
      assert(this->isInitialized());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() == this->outRows());
      // set sizes : Components/Transform/src/Fft/Backend/Fftw/ChebyshevEnergy.cpp
      std::cerr << "in =  \n";
      std::cerr << in;
      std::cerr << "\n";

      auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
      auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
      auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);
      this->applyPreOperator(tmpIn, in, true);
      this->mBackend.applyFft(tmpOut, tmpIn); // calculates toroidal/poloidal scalars in physical space (with pre-operator applied earlier). Puts it in tmpOut
      this->mBackend.square(tmpSquare, tmpOut, true); // squares it, put the result in tmpSquare
      this->applyPreOperator(tmpIn, in, false);
      this->mBackend.applyFft(tmpOut, tmpIn);
      this->mBackend.square(tmpSquare, tmpOut, false);
      this->mBackend.applyFwdFft(tmpOut, tmpSquare); // calculates Chebyshev coefficients of tmpSquare
      this->applyPostOperator(rOut, tmpOut); //this does the correct energy integral
      // For Energy.cpp, this integration is found in -> Fftw/ChebyshevEnergy.cpp::output(Matrix& rOut, const Matrix& tmp) 
      // which calculates the correct integral (NOT an FFT), using \int T_n dx (see for example Wikipedia)
   }

   // anelastic version:
   void IChebyshevEnergy::transform(Matrix& rOut, const MatrixZ& in, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
   {
      // for dim1D=8, dim2D=dim3D=4
      // blockSize = 15
      // rOut = 15 x 1
      // in = 26 x 15 (I think both anelastic and boussinesq)
      assert(this->isInitialized());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() == this->outRows());

      auto eGrid = this->mBackend.getEGrid();
      
      auto rho = pF->evaluate(eGrid,0,0);

      auto& tmpIn = this->mBackend.getStorage(StorageKind::in); // Boussinesq size : m_rows = 52, m_cols = 15; Anelastic: m_row = 84
      auto& tmpOut = this->mBackend.getStorage(StorageKind::out); // Boussinesq size : m_rows = 52, m_cols = 15; Anelastic: m_row = 84 
      auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid); // Boussinesq size : m_rows = 52, m_cols = 15; Anelastic: m_row = 84 
      // BUT ENERGYY2 is not calling the initBackend, uses Boussinesq size... (or does it.....)
      // SOMETHING IS OFF WITH THE SIZES: ANELASTIC = BOUSSINESQ....

      std::cerr << "in =  \n";
      std::cerr << in;
      std::cerr << "\n";

      // Set small values to zero
      //MatrixZ cleanedIn = in;
      //const MHDFloat threshold = 1e-100;
      //cleanedIn = (cleanedIn.array().abs() < threshold).select(0.0, cleanedIn);

      // Set small values to zero - more aggressive cleaning
      /*
      MatrixZ cleanedIn = in;
      const MHDFloat threshold = 1e-50;  // Much more aggressive threshold

      // Clean both real and imaginary parts separately
      for (int i = 0; i < cleanedIn.rows(); ++i) {
         for (int j = 0; j < cleanedIn.cols(); ++j) {
            auto& val = cleanedIn(i, j);
            if (std::abs(val.real()) < threshold) {
                  val = std::complex<MHDFloat>(0.0, val.imag());
            }
            if (std::abs(val.imag()) < threshold) {
                  val = std::complex<MHDFloat>(val.real(), 0.0);
            }
            if (std::abs(val) < threshold) {
                  val = std::complex<MHDFloat>(0.0, 0.0);
            }
         }
      }
      */
      //std::cerr << "in =  \n";
      //std::cerr << cleanedIn;
      //std::cerr << "\n";

      // in  contains e-310 numbers, not a problem for Boussinesq
      this->applyPreOperator(tmpIn, in, true); // tmpIn contains zeros and 1 in Boussinesq
                                               // Anelastic: contains e-310 and e+194  numbers !!!


      std::cerr << "tmpIn =  \n";
      std::cerr << tmpIn;
      std::cerr << "\n";



      this->mBackend.applyFft(tmpOut, tmpIn); // Boussinesq gives second col = 1
                                              // anelasitc gives second col = 1, bit also some e+175

      std::cerr << "tmpOut =  \n"; 
      std::cerr << tmpOut;
      std::cerr << "\n";


      this->mBackend.square(tmpSquare, tmpOut, true); // therefore gives inf

      std::cerr << "tmpSquare =  \n";
      std::cerr << tmpSquare;
      std::cerr << "\n";


      this->applyPreOperator(tmpIn, in, false);

      std::cerr << "tmpIn =  \n";
      std::cerr << tmpIn;
      std::cerr << "\n";

      this->mBackend.applyFft(tmpOut, tmpIn); // give too big numbers
      
      std::cerr << "tmpOut =  \n";
      std::cerr << tmpOut;
      std::cerr << "\n";

      this->mBackend.square(tmpSquare, tmpOut, false); // therefore gives inf

      std::cerr << "tmpSquare =  \n";
      std::cerr << tmpSquare;
      std::cerr << "\n";


      tmpSquare = tmpSquare.array().colwise() / rho.array(); // divides energy by rho
      this->mBackend.applyFwdFft(tmpOut, tmpSquare);
      this->applyPostOperator(rOut, tmpOut);

      std::cerr << "tmpOut =  \n";
      std::cerr << tmpOut;
      std::cerr << "\n";


   }


   void IChebyshevEnergy::transform(Matrix& rOut, const Matrix& in) const
   {
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
   }

   void IChebyshevEnergy::transform(Matrix& rOut, const Matrix& in, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
   {
      assert(this->isInitialized());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() == this->outRows());

      auto eGrid = this->mBackend.getEGrid();
      
      auto rho = pF->evaluate(eGrid,0,0);

      auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
      auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
      auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);
      this->applyPreOperator(tmpIn, in);
      this->mBackend.applyFft(tmpOut, tmpIn);
      this->mBackend.square(tmpSquare, tmpOut, true);
      tmpSquare = tmpSquare.array().colwise() / rho.array(); // divides energy by rho
      this->mBackend.applyFwdFft(tmpOut, tmpSquare);
      this->applyPostOperator(rOut, tmpOut);
   }

   void IChebyshevEnergy::transform(MatrixZ&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Chebyshev FFT energy reductor");
   }

   void IChebyshevEnergy::transform(MatrixZ&, const Matrix&) const
   {
      throw std::logic_error("Data is not compatible with Chebyshev FFT energy reductor");
   }

   MHDFloat IChebyshevEnergy::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   int IChebyshevEnergy::outRows() const
   {
      return this->mspSetup->blockSize();
   }

   int IChebyshevEnergy::outCols() const
   {
      return 1;
   }

}
}
}
}
}
}
