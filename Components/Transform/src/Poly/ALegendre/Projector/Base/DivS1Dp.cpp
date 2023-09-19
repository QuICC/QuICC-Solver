/**
 * @file DivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre 1/S D_phi projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/DivS1Dp.hpp"
#include "Types/Constants.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   void DivS1Dp<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      MHDComplex c(0.0,static_cast<MHDFloat>(this->mspSetup->slow(i)));
      DivS1<base_t>::applyOperator(rOut, i, c*in);
   }

}
}
}
}
}
