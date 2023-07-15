/**
 * @file Ll2.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1)^2 P integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/Ll2.hpp"
#include "QuICC/Polynomial/ALegendre/Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   void Ll2<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      P<base_t>::applyOperator(rOut, i, in);
      rOut = this->mLl2.bottomRows(rOut.rows()).asDiagonal()*rOut;
   }

   void Ll2<base_t>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLl2 = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl2 = (this->mLl2.array()*(this->mLl2.array() + 1.0)).pow(2);
   }

}
}
}
}
}
