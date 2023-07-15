/**
 * @file LlDivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1)/Sin D_phi projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/LlDivS1Dp.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   void LlDivS1Dp<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      DivS1::applyOperator(rOut, i, (Math::cI*static_cast<MHDFloat>(this->mspSetup->slow(i))*this->mLl.bottomRows(in.rows())).asDiagonal()*in);
   }

   void LlDivS1Dp<base_t>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl = this->mLl.array()*(this->mLl.array() + 1.0);
   }

}
}
}
}
}
