/**
 * @file Power.cpp
 * @brief Source of the implementation of the Bessel power spectrum operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/Power.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   Power<base_t>::Power()
      : Bessel::Reductor::IBesselPower(1)
   {
      this->setProfileTag();
   }

   void Power<base_t>::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
