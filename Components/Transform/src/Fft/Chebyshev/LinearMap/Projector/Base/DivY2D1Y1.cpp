/**
 * @file DivY2D1Y1.cpp
 * @brief Source of the implementation of the Chebyshev 1/Y^2 D Y projector, with linear map y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y2.hpp"
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/Base/DivY2D1Y1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {


   void DivY2D1Y1<base_t>::initOperator() const
   {
      // Check for division by 0!
      assert(this->mspSetup->lower() > 0.0 || this->mspSetup->upper() < 0.0);

      // For (1/y^2)D1Y1, we solve Y1 f = I1 g and then apply the 1/y^2 scaling
      // First set up I1 operator
      ::QuICC::SparseSM::Chebyshev::LinearMap::I1 op(this->mspSetup->specSize()+2,
                                                     this->mspSetup->specSize()+2, 
                                                     this->mspSetup->lower(), 
                                                     this->mspSetup->upper());

      this->mBackend.solver().setOperator(op.mat(), 1, 1);

      // Set up Y1 operator on left side
      ::QuICC::SparseSM::Chebyshev::LinearMap::Y1 opY1(this->mspSetup->specSize()+1, // Maps N to N+1
                                                       this->mspSetup->specSize(), 
                                                       this->mspSetup->lower(), 
                                                       this->mspSetup->upper());

      this->mBackend.solver().setSpectralOperator(opY1.mat(), 1);

      // Set up Chebyshev quadrature for physical space operations
      Internal::Array igrid, iweights;
      Polynomial::Quadrature::ChebyshevRule quad;
      quad.computeQuadrature(igrid, iweights, this->mspSetup->fwdSize(), this->mspSetup->lower(), this->mspSetup->upper());
      // Set up 1/y^2 scaling for physical space
      this->mBackend.setScaler(igrid.array().pow(-2).cast<MHDFloat>().matrix());
   }

   void DivY2D1Y1<base_t>::initBackend() const
   {
      // Call parent initializer
      ILinearMapProjector::initBackend();

      // Initialize the solver with 1 extra mode for first derivative
      this->mBackend.addSolver(1);
   }

   void DivY2D1Y1<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const Matrix>& in) const
   {
      this->mBackend.input(tmp, in);
      // Apply spectral operator (Y1)
      auto specOp = this->mBackend.solver().getSpectralOperator();
      tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
      // Solve system Y1 f = I1 g
      this->mBackend.getSolution(tmp, 1, 1);
   }

   void DivY2D1Y1<base_t>::applyPostOperator(Eigen::Ref<Matrix> rOut) const
   {
       // Apply 1/y^2 scaling in physical space
      this->mBackend.outputScale(rOut);
   }

   // Complex version of pre-operator
   void DivY2D1Y1<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const MatrixZ>& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
      auto specOp = this->mBackend.solver().getSpectralOperator();
      tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
      this->mBackend.getSolution(tmp, 1, 1);
   }

   // Complex version of post-operator
   void DivY2D1Y1<base_t>::applyPostOperator(Eigen::Ref<MatrixZ> rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mBackend.outputScale(rOut, tmp, useReal);
   }

}
}
}
}
}
}
