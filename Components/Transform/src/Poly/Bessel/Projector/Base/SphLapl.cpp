/**
 * @file SphLapl.cpp
 * @brief Source of the implementation of the Bessel spherical laplacian projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Projector/Base/SphLapl.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Projector {

   SphLapl<base_t>::SphLapl()
   {
      this->setProfileTag();
   }

   void SphLapl<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
