/**
 * @file I2D1_I2.cpp
 * @brief Source of the implementation of the Chebyshev I^2 of D integrator, but 0 mode is I^2 of P integrator, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2D1_I2.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   I2D1_I2::I2D1_I2()
   {
   }

   I2D1_I2::~I2D1_I2()
   {
   }

   void I2D1_I2::initOperator() const
   {
      int size = this->mspSetup->specSize() + std::min(1, this->mspSetup->padSize());
      ::QuICC::SparseSM::Chebyshev::LinearMap::I2D1 op(size, size, this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setSpectralOperator(op.mat().topRows(this->mspSetup->specSize()));

      if(this->mspSetup->slowSize() > 0)
      {
         size = this->mspSetup->specSize() + std::min(2, this->mspSetup->padSize());
         ::QuICC::SparseSM::Chebyshev::LinearMap::I2 meanOp(size, size, this->mspSetup->lower(), this->mspSetup->upper());
         this->mBackend.setMeanOperator(meanOp.mat().topRows(this->mspSetup->specSize()));
      }
   }

   void I2D1_I2::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.outputSpectral(rOut);
   }

   void I2D1_I2::applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
   }

   void I2D1_I2::applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mBackend.outputSpectral(rOut, tmp, useReal);
   }

}
}
}
}
}
}
