/**
 * @file Ll.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1) P projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/Ll.hpp"
#include "QuICC/Polynomial/ALegendre/Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   void Ll<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      P<base_t>::applyOperator(rOut, i, this->mLl.bottomRows(in.rows()).asDiagonal()*in);
   }

   void Ll<base_t>::initSpecial() const
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
