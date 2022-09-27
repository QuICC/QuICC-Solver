/**
 * @file DivY1.cpp
 * @brief Source of the implementation of the Chebyshev 1/Y projector, with linear map y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1.hpp"

// Project includes
//
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   DivY1::DivY1()
   {
   }

   DivY1::~DivY1()
   {
   }

   void DivY1::initOperator() const
   {
      // Check for division by 0!
      assert(this->mspSetup->lower() > 0.0 || this->mspSetup->upper() < 0.0);

      internal::Array igrid, iweights;
      Polynomial::Quadrature::ChebyshevRule quad;
      quad.computeQuadrature(igrid, iweights, this->mspSetup->fwdSize(), this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setScaler(igrid.array().pow(-1).cast<MHDFloat>().matrix());
   }

   void DivY1::applyPreOperator(Matrix& tmp, const Matrix& in) const
   {
      this->mBackend.input(tmp, in);
   }

   void DivY1::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.outputScale(rOut);
   }

   void DivY1::applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
   }

   void DivY1::applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mBackend.outputScale(rOut, tmp, useReal);
   }

}
}
}
}
}
}
