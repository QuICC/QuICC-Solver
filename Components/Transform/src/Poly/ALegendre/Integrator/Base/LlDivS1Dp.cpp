/**
 * @file LlDivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1) 1/Sin D_phi integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/LlDivS1Dp.hpp"
#include "Types/Math.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   void LlDivS1Dp<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      DivS1<base_t>::applyOperator(rOut, i, in);
      rOut = (static_cast<MHDFloat>(-this->mspSetup->slow(i))*Math::cI*this->mLl.bottomRows(rOut.rows())).asDiagonal()*rOut;
   }

   void LlDivS1Dp<base_t>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl = (this->mLl.array()*(this->mLl.array() + 1.0));
   }

}
}
}
}
}
