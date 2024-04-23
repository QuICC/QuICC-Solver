/**
 * @file P_Zero.cpp
 * @brief Source of the implementation of the Bessel P_Zero integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Integrator/Base/P_Zero.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Integrator {

   P_Zero<base_t>::P_Zero()
      : IBesselIntegrator()
   {
      this->setProfileTag();
   }

   void P_Zero<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
