/**
 * @file DivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre 1/Sin D_phi integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/DivS1Dp.hpp"
#include "Types/Constants.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   void DivS1Dp<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      DivS1<base_t>::applyOperator(rOut, i, in);
      rOut = static_cast<MHDFloat>(-this->mspSetup->slow(i))*Math::cI*rOut;
   }

}
}
}
}
}
