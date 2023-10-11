/**
 * @file SphRadLapl.cpp
 * @brief Source of the implementation of the Chebyshev radial part of spherical laplacian 1/Y^2 D Y^2 projector, with linear map y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/SphRadLapl.hpp"

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y2.hpp"
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   SphRadLapl::SphRadLapl()
   {
   }

   SphRadLapl::~SphRadLapl()
   {
   }

   void SphRadLapl::initOperator() const
   {
      // Check for division by 0!
      assert(this->mspSetup->lower() > 0.0 || this->mspSetup->upper() < 0.0);

      ::QuICC::SparseSM::Chebyshev::LinearMap::I1 spasmI1(this->mspSetup->specSize()+3,this->mspSetup->specSize()+3, this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.solver().setOperator(spasmI1.mat(), 2, 2);

      ::QuICC::SparseSM::Chebyshev::LinearMap::Y2 spasmY2(this->mspSetup->specSize()+2,this->mspSetup->specSize()+2, this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.solver().setSpectralOperator(spasmY2.mat(), 1);

      Internal::Array igrid, iweights;
      Polynomial::Quadrature::ChebyshevRule quad;
      quad.computeQuadrature(igrid, iweights, this->mspSetup->fwdSize(), this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setScaler(igrid.array().pow(-2).cast<MHDFloat>().matrix());
   }

   void SphRadLapl::initBackend() const
   {
      // Call parent initializer
      IChebyshevProjector::initBackend();

      // Initialize the solver
      this->mBackend.addSolver(2);
   }

   void SphRadLapl::applyPreOperator(Matrix& tmp, const Matrix& in) const
   {
      this->mBackend.input(tmp, in, 1);
      this->mBackend.getSolution(tmp, 3, -1, true);
      auto specOp = this->mBackend.solver().getSpectralOperator();
      tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
      this->mBackend.getSolution(tmp, 1, 2);
   }

   void SphRadLapl::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.outputScale(rOut);
   }

   void SphRadLapl::applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, 1, useReal);
      this->mBackend.getSolution(tmp, 3, -1, true);
      auto specOp = this->mBackend.solver().getSpectralOperator();
      tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
      this->mBackend.getSolution(tmp, 1, 2);
   }

   void SphRadLapl::applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mBackend.outputScale(rOut, tmp, useReal);
   }

}
}
}
}
}
}
