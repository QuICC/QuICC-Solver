/**
 * @file LlD1.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1) D projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/LlD1.hpp"
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   void LlD1<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      D1<base_t>::applyOperator(rOut, i, this->mLl.bottomRows(in.rows()).asDiagonal()*in);
   }

   void LlD1<base_t>::initSpecial() const
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
