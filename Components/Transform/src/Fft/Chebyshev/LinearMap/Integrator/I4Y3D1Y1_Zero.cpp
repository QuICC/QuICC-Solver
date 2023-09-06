/**
 * @file I4Y3D1Y1_Zero.cpp
 * @brief Source of the implementation of the  Chebyshev I^4 Y^4 of 1/Y D Y integrator, but 0 mode is zeroed, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4Y3D1Y1_Zero.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y3D1Y1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   I4Y3D1Y1_Zero::I4Y3D1Y1_Zero()
   {
   }

   I4Y3D1Y1_Zero::~I4Y3D1Y1_Zero()
   {
   }

   void I4Y3D1Y1_Zero::initOperator() const
   {
      int size = this->mspSetup->specSize() + std::min(7, this->mspSetup->padSize());
      ::QuICC::SparseSM::Chebyshev::LinearMap::I4Y3D1Y1 op(size, size, this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setSpectralOperator(op.mat().topRows(this->mspSetup->specSize()));

      if(this->mspSetup->slowSize() > 0)
      {
         this->mBackend.setMeanOperator(SparseMatrix(size, size).topRows(this->mspSetup->specSize()));
      }
   }

   void I4Y3D1Y1_Zero::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.outputSpectral(rOut);
   }

   void I4Y3D1Y1_Zero::applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
   }

   void I4Y3D1Y1_Zero::applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mBackend.outputSpectral(rOut, tmp, useReal);
   }

}
}
}
}
}
}
