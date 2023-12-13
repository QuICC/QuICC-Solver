/**
 * @file PDivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre 1/Sin D_phi parallel integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/DivS1Dp.hpp"
#include "Types/Math.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   void DivS1Dp<kokkos_t>::makeOperator(Matrix &op,
      const Internal::Array &igrid, const Internal::Array &iweights, const int i) const {
       DivS1<kokkos_t>::makeOperator(op, igrid, iweights, i);
       op = static_cast<MHDFloat>(-this->mspSetup->slow(i)) * op;
   }

   void DivS1Dp<kokkos_t>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      DivS1<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);

      DataType constant;
      constant.real() = Math::cI.real();
      constant.imag() = Math::cI.imag();

      constantMultiplyMatrix(constant, rOutView);
   }

}
}
}
}
}
