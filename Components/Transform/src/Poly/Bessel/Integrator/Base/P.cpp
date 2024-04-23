/**
 * @file P.cpp
 * @brief Source of the implementation of the Bessel P integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Integrator/Base/P.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Integrator {

   P<base_t>::P()
      : IBesselIntegrator()
   {
      this->setProfileTag();
   }

   void P<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
