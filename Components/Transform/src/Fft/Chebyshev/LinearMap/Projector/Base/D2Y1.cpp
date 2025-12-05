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

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/Base/D2Y1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   void D2Y1<base_t>::initOperator() const
   {
      // We want to find y from D2Y1<base_t> x = y
      // To solve this: I2D2 (Y1 x) = I2 y, where I2D2=identity with 2 zero rows.
      // Y1 x = I2 y can be inverted for y
      //
      // The multiplication by Y1 increases the expansion by 1, we need the differntial 
      // operator of size (N+1) x (N+1)
      // I2 has two zero rows we need to compute I2 of size (N+3) x (N+3) in order to extra
      // The bottom left sub matrix of size (N+1) x (N+1)
      ::QuICC::SparseSM::Chebyshev::LinearMap::I2 op(this->mspSetup->specSize()+3,this->mspSetup->specSize()+3, this->mspSetup->lower(), this->mspSetup->upper());

      // The last two args of setOperator are extraRows and extraColumns (see DifferentialSolver.cpp)
      // We extract the submatrix from bottom left corner of size N+1 x N+1
      this->mBackend.solver().setOperator(op.mat(), 1, 1);

      // Y1 increases the spectral expansion by 1 we need size (N+1) x N
      // We dropped the first two rows in I2, so we need to shift Y1 by 2
      ::QuICC::SparseSM::Chebyshev::LinearMap::Y1 opY1(this->mspSetup->specSize()+1,this->mspSetup->specSize(), this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.solver().setSpectralOperator(opY1.mat(), 2);
   }

   void D2Y1<base_t>::initBackend() const
   {
      // Call parent initializer
      ILinearMapProjector::initBackend();

      // Initialize the solver
      this->mBackend.addSolver(1); // I2 has size N+1
   }

   void D2Y1<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const Matrix>& in) const
   {
      this->mBackend.input(tmp, in);
      auto specOp = this->mBackend.solver().getSpectralOperator();
      tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
      this->mBackend.getSolution(tmp, 2, 1); // Last 2 coefficients are zero (D^2)
                                             // 1 extra modes
   }

   void D2Y1<base_t>::applyPostOperator(Eigen::Ref<Matrix>) const
   {
   }

   void D2Y1<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const MatrixZ>& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
      auto specOp = this->mBackend.solver().getSpectralOperator();
      tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
      this->mBackend.getSolution(tmp, 2, 1); // Last 2 coefficients are zero (D^2)
                                             // 1 extra mode
   }

   void D2Y1<base_t>::applyPostOperator(Eigen::Ref<MatrixZ> rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mBackend.output(rOut, tmp, useReal);
   }

}
}
}
}
}
}
