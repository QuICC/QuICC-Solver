/**
 * @file I4.cpp
 * @brief Source of the implementation of the Chebyshev I^4 of P integrator, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   I4::I4()
   {
   }

   I4::~I4()
   {
   }

   void I4::initOperator() const
   {
      int size = this->mspSetup->specSize() + std::min(4, this->mspSetup->padSize());
      ::QuICC::SparseSM::Chebyshev::LinearMap::I4 op(size, size, this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setSpectralOperator(op.mat().topRows(this->mspSetup->specSize()));
   }

   void I4::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.outputSpectral(rOut);
   }

   void I4::applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
   }

   void I4::applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mBackend.outputSpectral(rOut, tmp, useReal);
   }

}
}
}
}
}
}
