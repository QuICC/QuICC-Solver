/**
 * @file D4.cpp
 * @brief Source of the implementation of the Chebyshev D^4 projector, with linear map y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D4.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   D4::D4()
   {
   }

   D4::~D4()
   {
   }

   void D4::initOperator() const
   {
      ::QuICC::SparseSM::Chebyshev::LinearMap::I4 op(this->mspSetup->specSize()+4,this->mspSetup->specSize()+4, this->mspSetup->lower(), this->mspSetup->upper());

      this->mBackend.solver().setOperator(op.mat());
   }

   void D4::initBackend() const
   {
      // Call parent initializer
      IChebyshevProjector::initBackend();

      // Initialize the solver
      this->mBackend.addSolver();
   }

   void D4::applyPreOperator(Matrix& rOut, const Matrix& in) const
   {
      this->mBackend.solver().input(in, 4);

      this->mBackend.getSolution(4);

      this->mBackend.output(rOut);
   }

   void D4::applyPostOperator(Matrix&) const
   {
   }

   void D4::applyPreOperator(const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.solver().input(in, 4, useReal);

      this->mBackend.getSolution(4);

      this->mBackend.io();
   }

   void D4::applyPostOperator(MatrixZ& rOut, const bool useReal) const
   {
      this->mBackend.output(rOut, useReal);
   }

}
}
}
}
}
}
