/**
 * @file DivY2.cpp
 * @brief Source of the implementation of the Chebyshev 1/R^2 projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY2.hpp"

// Project includes
//
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   DivY2::DivY2()
   {
   }

   DivY2::~DivY2()
   {
   }

   void DivY2::initOperator() const
   {
      // Check for division by 0!
      assert(this->mspSetup->lower() > 0.0 || this->mspSetup->upper() < 0.0);

      internal::Array igrid, iweights;
      Polynomial::Quadrature::ChebyshevRule quad;
      quad.computeQuadrature(igrid, iweights, this->mspSetup->fwdSize(), this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setScaler(igrid.array().pow(-2).cast<MHDFloat>().matrix());
   }

   void DivY2::applyPreOperator(Matrix& rOut, const Matrix& in) const
   {
      this->mBackend.input(in, true);

      this->mBackend.output(rOut);
   }

   void DivY2::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.outputScale(rOut);
   }

   void DivY2::applyPreOperator(const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(in, useReal, true);

      this->mBackend.io();
   }

   void DivY2::applyPostOperator(MatrixZ& rOut, const bool useReal) const
   {
      this->mBackend.outputScale(rOut, useReal);
   }

}
}
}
}
}
}
