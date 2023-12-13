/**
 * @file DivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre 1/S D_phi projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/DivS1Dp.hpp"
#include "Types/Math.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   void DivS1Dp<kokkos_t>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      constantMultiplyMatrix<1>(this->mspSetup, scan, inView);
      DivS1<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);
   }

}
}
}
}
}

