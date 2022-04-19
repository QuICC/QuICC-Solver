/**
 * @file I2_I2D1.cpp
 * @brief Source of the implementation of the Chebyshev I^2 of P integrator, but 0 mode is I^2 of D integrator, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2_I2D1.hpp"

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

   I2_I2D1::I2_I2D1()
   {
   }

   I2_I2D1::~I2_I2D1()
   {
   }

   void I2_I2D1::initOperator() const
   {
      int size = this->mspSetup->specSize() + std::min(2, this->mspSetup->padSize());
      ::QuICC::SparseSM::Chebyshev::LinearMap::I2 op(size, size, this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setSpectralOperator(op.mat().topRows(this->mspSetup->specSize()));

      if(this->mspSetup->slowSize() > 0)
      {
         size = this->mspSetup->specSize() + std::min(1, this->mspSetup->padSize());
         ::QuICC::SparseSM::Chebyshev::LinearMap::I2D1 meanOp(size, size, this->mspSetup->lower(), this->mspSetup->upper());
         this->mBackend.setMeanOperator(op.mat().topRows(this->mspSetup->specSize()));
      }
   }

   void I2_I2D1::applyPreOperator(Matrix& rOut, const Matrix& in) const
   {
      this->mBackend.io(rOut, in);
   }

   void I2_I2D1::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.outputSpectral(rOut);
   }

   void I2_I2D1::applyPreOperator(const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(in, useReal);
   }

   void I2_I2D1::applyPostOperator(MatrixZ& rOut, const bool useReal) const
   {
      this->mBackend.outputSpectral(rOut, useReal);
   }

}
}
}
}
}
}
