/**
 * @file P.cpp
 * @brief Source of the implementation of the Bessel P projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Projector/Base/P.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Projector {

   P<base_t>::P()
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
