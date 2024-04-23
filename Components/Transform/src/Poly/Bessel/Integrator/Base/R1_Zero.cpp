/**
 * @file R1_Zero.cpp
 * @brief Source of the implementation of the Bessel R1_Zero integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Integrator/Base/R1_Zero.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Integrator {

   R1_Zero<base_t>::R1_Zero()
      : IBesselIntegrator()
   {
      this->setProfileTag();
   }

   void R1_Zero<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
