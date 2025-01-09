/**
 * @file D2Y1.cpp
 * @brief Source of the implementation of the Chebyshev D^2 Y projector, with linear map y = ax + b
 * 
 * Modified from the D1Y1 projector
 * 
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D2Y1.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   D2Y1::D2Y1()
   {
   }

   D2Y1::~D2Y1()
   {
   }

   void D2Y1::initOperator() const
   {
      // We want to find y from D2Y2 x = y
      // To solve this: I2D2 (Y2 x) = I2 y, where I2D2=identity.
      // Y2 x = I2 y can be inverted for y
      //
      // In the argument of the specSize function,  specSize()+3 represents the bandwidth of the matrix
      // For I2 it should be specSize()+2, but the Y1 makes it a specSize()+3
      ::QuICC::SparseSM::Chebyshev::LinearMap::I2 op(this->mspSetup->specSize()+3,this->mspSetup->specSize()+3, this->mspSetup->lower(), this->mspSetup->upper());

      // The last two args of setOperator are extraRows and extraColumns (see DifferentialSolver.cpp)
      // Here they are set to 2 because of the order of I2.
      this->mBackend.solver().setOperator(op.mat(), 2, 2);

      ::QuICC::SparseSM::Chebyshev::LinearMap::Y1 opY1(this->mspSetup->specSize()+1,this->mspSetup->specSize(), this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.solver().setSpectralOperator(opY1.mat(), 1);
   }

   void D2Y1::initBackend() const
   {
      // Call parent initializer
      IChebyshevProjector::initBackend();

      // Initialize the solver
      this->mBackend.addSolver(1); // *** I think this is 1 (extra spectral degrees)
   }

   void D2Y1::applyPreOperator(Matrix& tmp, const Matrix& in) const
   {
      this->mBackend.input(tmp, in);
      auto specOp = this->mBackend.solver().getSpectralOperator();
      tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
      this->mBackend.getSolution(tmp, 1, 1); // **** Are these 2,2?
   }

   void D2Y1::applyPostOperator(Matrix&) const
   {
   }

   void D2Y1::applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
      auto specOp = this->mBackend.solver().getSpectralOperator();
      tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
      this->mBackend.getSolution(tmp, 1, 1);
   }

   void D2Y1::applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mBackend.output(rOut, tmp, useReal);
   }

}
}
}
}
}
}
