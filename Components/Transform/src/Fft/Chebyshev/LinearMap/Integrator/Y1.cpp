/**
 * @file Y.cpp
 * @brief Source of the implementation of the  Chebyshev Y integrator, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/Y1.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   Y1::Y1()
   {
   }

   Y1::~Y1()
   {
   }

   void Y1::initOperator() const
   {
      int size = this->mspSetup->specSize() + std::min(1, this->mspSetup->padSize());
      ::QuICC::SparseSM::Chebyshev::LinearMap::Y1 op(size, size, this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setSpectralOperator(op.mat().topRows(this->mspSetup->specSize()));
   }

   void Y1::applyPreOperator(Matrix& rOut, const Matrix& in) const
   {
      this->mBackend.io(rOut, in);
   }

   void Y1::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.outputSpectral(rOut);
   }

   void Y1::applyPreOperator(const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(in, useReal);
   }

   void Y1::applyPostOperator(MatrixZ& rOut, const bool useReal) const
   {
      this->mBackend.outputSpectral(rOut, useReal);
   }

}
}
}
}
}
}
