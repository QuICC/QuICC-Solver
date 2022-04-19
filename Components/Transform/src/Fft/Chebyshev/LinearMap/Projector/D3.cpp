/**
 * @file D3.cpp
 * @brief Source of the implementation of the Chebyshev D^3 projector, with linear map y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D3.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I3.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   D3::D3()
   {
   }

   D3::~D3()
   {
   }

   void D3::initOperator() const
   {
      ::QuICC::SparseSM::Chebyshev::LinearMap::I3 op(this->mspSetup->specSize()+3,this->mspSetup->specSize()+3, this->mspSetup->lower(), this->mspSetup->upper());

      this->mBackend.solver().setOperator(op.mat());
   }

   void D3::initBackend() const
   {
      // Call parent initializer
      IChebyshevProjector::initBackend();

      // Initialize the solver
      this->mBackend.addSolver();
   }

   void D3::applyPreOperator(Matrix& rOut, const Matrix& in) const
   {
      this->mBackend.solver().input(in, 3);

      this->mBackend.getSolution(3);

      this->mBackend.output(rOut);
   }

   void D3::applyPostOperator(Matrix&) const
   {
   }

   void D3::applyPreOperator(const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.solver().input(in, 3, useReal);

      this->mBackend.getSolution(3);

      this->mBackend.io();
   }

   void D3::applyPostOperator(MatrixZ& rOut, const bool useReal) const
   {
      this->mBackend.output(rOut, useReal);
   }

}
}
}
}
}
}
