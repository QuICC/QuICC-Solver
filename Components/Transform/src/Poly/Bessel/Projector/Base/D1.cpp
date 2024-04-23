/**
 * @file D1.cpp
 * @brief Source of the implementation of the Bessel D1 projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Projector/Base/D1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Projector {

   D1<base_t>::D1()
   {
      this->setProfileTag();
   }

   void D1<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
