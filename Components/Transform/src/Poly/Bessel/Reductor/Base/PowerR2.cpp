/**
 * @file PowerR2.cpp
 * @brief Source of the implementation of the Bessel R^2 power spectrum operator
 */

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/PowerR2.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   PowerR2<base_t>::PowerR2()
      : Bessel::Reductor::IBesselPower(1)
   {
      this->setProfileTag();
   }

   void PowerR2<base_t>::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
