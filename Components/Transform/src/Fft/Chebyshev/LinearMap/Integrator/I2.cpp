/**
 * @file I2.cpp
 * @brief Source of the implementation of the Chebyshev I2 of P integrator, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   I2::I2()
   {
   }

   I2::~I2()
   {
   }

   void I2::initOperator() const
   {
      int size = this->mspSetup->specSize() + std::min(2, this->mspSetup->padSize());
      ::QuICC::SparseSM::Chebyshev::LinearMap::I2 op(size, size, this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setSpectralOperator(op.mat().topRows(this->mspSetup->specSize()));
   }

   void I2::applyPreOperator(Matrix& rOut, const Matrix& in) const
   {
      this->mBackend.io(rOut, in);
   }

   void I2::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.outputSpectral(rOut);
   }

   void I2::applyPreOperator(const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(in, useReal);
   }

   void I2::applyPostOperator(MatrixZ& rOut, const bool useReal) const
   {
      this->mBackend.outputSpectral(rOut, useReal);
   }

}
}
}
}
}
}
