/**
 * @file DivLl.cpp
 * @brief Source of the implementation of the associated Legendre 1/l(l+1) P integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/DivLl.hpp"
#include "QuICC/Polynomial/ALegendre/Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   void DivLl<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      P<base_t>::applyOperator(rOut, i, in);
      rOut = this->mDivLl.bottomRows(rOut.rows()).asDiagonal()*rOut;
   }

   void DivLl<base_t>::initSpecial() const
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
