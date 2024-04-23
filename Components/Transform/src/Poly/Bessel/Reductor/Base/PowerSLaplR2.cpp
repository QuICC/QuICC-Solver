/**
 * @file PowerSLaplR2.cpp
 * @brief Source of the implementation of the Bessel Spherical Laplacian R^2 power spectrum operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/PowerSLaplR2.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   PowerSLaplR2<base_t>::PowerSLaplR2()
      : Bessel::Reductor::IBesselPower(1)
   {
      this->setProfileTag();
   }

   void PowerSLaplR2<base_t>::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
