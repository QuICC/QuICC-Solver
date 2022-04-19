/**
 * @file D1.cpp
 * @brief Source of the implementation of the Chebyshev D projector, with linear map y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   D1::D1()
   {
   }

   D1::~D1()
   {
   }

   void D1::initOperator() const
   {
      ::QuICC::SparseSM::Chebyshev::LinearMap::I1 op(this->mspSetup->specSize()+1,this->mspSetup->specSize()+1, this->mspSetup->lower(), this->mspSetup->upper());

      this->mBackend.solver().setOperator(op.mat());
   }

   void D1::initBackend() const
   {
      // Call parent initializer
      IChebyshevProjector::initBackend();

      // Initialize the solver
      this->mBackend.addSolver();
   }

   void D1::applyPreOperator(Matrix& rOut, const Matrix& in) const
   {
      this->mBackend.solver().input(in, 1);

      this->mBackend.getSolution(1);

      this->mBackend.output(rOut);
   }

   void D1::applyPostOperator(Matrix&) const
   {
   }

   void D1::applyPreOperator(const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.solver().input(in, 1, useReal);

      this->mBackend.getSolution(1);

      this->mBackend.io();
   }

   void D1::applyPostOperator(MatrixZ& rOut, const bool useReal) const
   {
      this->mBackend.output(rOut, useReal);
   }

}
}
}
}
}
}
