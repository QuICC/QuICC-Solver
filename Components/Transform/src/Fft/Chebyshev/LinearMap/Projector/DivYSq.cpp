/**
 * @file DivYSq.cpp
 * @brief Source of the implementation of the Chebyshev 1/R^2 projector
 * equivalent to DivY2
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivYSq.hpp"

// Project includes
//
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   DivYSq::DivYSq()
   {
   }

   DivYSq::~DivYSq()
   {
   }

   void DivYSq::initOperator() const
   {
      // Check for division by 0!
      assert(this->mspSetup->lower() > 0.0 || this->mspSetup->upper() < 0.0);

      Internal::Array igrid, iweights;
      Polynomial::Quadrature::ChebyshevRule quad;
      quad.computeQuadrature(igrid, iweights, this->mspSetup->fwdSize(), this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setScaler(igrid.array().pow(-2).cast<MHDFloat>().matrix());
   }

   void DivYSq::applyPreOperator(Matrix& tmp, const Matrix& in) const
   {
      this->mBackend.input(tmp, in);
   }

   void DivYSq::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.outputScale(rOut);
   }

   void DivYSq::applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
   }

   void DivYSq::applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mBackend.outputScale(rOut, tmp, useReal);
   }

}
}
}
}
}
}
