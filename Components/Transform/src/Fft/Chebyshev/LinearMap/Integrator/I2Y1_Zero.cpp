/**
 * @file I2Y1_Zero.cpp
 * @brief Source of the implementation of the  Chebyshev I^2 Y^2 of 1/Y integrator, but 0 mode is zeroed, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2Y1_Zero.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   I2Y1_Zero::I2Y1_Zero()
   {
   }

   I2Y1_Zero::~I2Y1_Zero()
   {
   }

   void I2Y1_Zero::initOperator() const
   {
      int size = this->mspSetup->specSize() + std::min(2, this->mspSetup->padSize());
      ::QuICC::SparseSM::Chebyshev::LinearMap::I2Y1 op(size, size, this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setSpectralOperator(op.mat().topRows(this->mspSetup->specSize()));

      if(this->mspSetup->slowSize() > 0)
      {
         this->mBackend.setMeanOperator(SparseMatrix(size, size).topRows(this->mspSetup->specSize()));
      }
   }

   void I2Y1_Zero::applyPreOperator(Matrix& rOut, const Matrix& in) const
   {
      this->mBackend.io(rOut, in);
   }

   void I2Y1_Zero::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.outputSpectral(rOut);
   }

   void I2Y1_Zero::applyPreOperator(const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(in, useReal);
   }

   void I2Y1_Zero::applyPostOperator(MatrixZ& rOut, const bool useReal) const
   {
      this->mBackend.outputSpectral(rOut, useReal);
   }

}
}
}
}
}
}
