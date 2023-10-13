/**
 * @file DivLlDivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre 1/l(l+1) 1/Sin D_phi integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/DivLlDivS1Dp.hpp"
#include "Types/Math.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   void DivLlDivS1Dp<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      DivS1::applyOperator(rOut, i, in);
      MHDComplex c(0.0,-static_cast<MHDFloat>(this->mspSetup->slow(i)));
      rOut = (c*this->mDivLl.bottomRows(rOut.rows())).asDiagonal()*rOut;
   }

   void DivLlDivS1Dp<base_t>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mDivLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mDivLl = (this->mDivLl.array()*(this->mDivLl.array() + 1.0)).pow(-1);
      this->mDivLl(0) = 0.0;
   }

}
}
}
}
}
