/**
 * @file EnergyY2.cpp
 * @brief Source of the implementation of the Chebyshev energy Y^2 reductor, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/EnergyY2.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   EnergyY2::EnergyY2()
   {
   }

   EnergyY2::~EnergyY2()
   {
   }

   void EnergyY2::initOperator() const
   {
      int size = 2*this->mspSetup->specSize() + std::min(2, 2*this->mspSetup->padSize());
      ::QuICC::SparseSM::Chebyshev::LinearMap::Y2 op(size, size, this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setSpectralOperator(op.mat().leftCols(2*this->mspSetup->specSize()));
   }

   void EnergyY2::applyPreOperator(const Matrix& in) const
   {
      this->mBackend.input(in, true);

      this->mBackend.io();
   }

   void EnergyY2::applyPostOperator(Matrix& rOut) const
   {
      assert(rOut.cols() == 1);
      this->mBackend.outputSpectral(rOut);
   }

   void EnergyY2::applyPreOperator(const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(in, useReal, true);

      this->mBackend.io();
   }

}
}
}
}
}
}
