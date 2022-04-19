/**
 * @file D2.cpp
 * @brief Source of the implementation of the Chebyshev D^2 projector, with linear map y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D2.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   D2::D2()
   {
   }

   D2::~D2()
   {
   }

   void D2::initOperator() const
   {
      ::QuICC::SparseSM::Chebyshev::LinearMap::I2 op(this->mspSetup->specSize()+2,this->mspSetup->specSize()+2, this->mspSetup->lower(), this->mspSetup->upper());

      this->mBackend.solver().setOperator(op.mat());
   }

   void D2::initBackend() const
   {
      // Call parent initializer
      IChebyshevProjector::initBackend();

      // Initialize the solver
      this->mBackend.addSolver();
   }

   void D2::applyPreOperator(Matrix& rOut, const Matrix& in) const
   {
      this->mBackend.solver().input(in, 2);

      this->mBackend.getSolution(2);

      this->mBackend.output(rOut);
   }

   void D2::applyPostOperator(Matrix&) const
   {
   }

   void D2::applyPreOperator(const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.solver().input(in, 2, useReal);

      this->mBackend.getSolution(2);

      this->mBackend.io();
   }

   void D2::applyPostOperator(MatrixZ& rOut, const bool useReal) const
   {
      this->mBackend.output(rOut, useReal);
   }

}
}
}
}
}
}
