/**
 * @file PowerD1R1.cpp
 * @brief Source of the implementation of the Bessel D R power spectrum operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/PowerD1R1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   PowerD1R1<base_t>::PowerD1R1()
      : Bessel::Reductor::IBesselPower(0)
   {
      this->setProfileTag();
   }

   void PowerD1R1<base_t>::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
