/**
 * @file DivR1D1R1_Zero.cpp
 * @brief Source of the implementation of the Bessel 1/R D R projector and zero for l = 0
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Projector/Base/DivR1D1R1_Zero.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Projector {

   DivR1D1R1_Zero<base_t>::DivR1D1R1_Zero()
   {
      this->setProfileTag();
   }

   void DivR1D1R1_Zero<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
