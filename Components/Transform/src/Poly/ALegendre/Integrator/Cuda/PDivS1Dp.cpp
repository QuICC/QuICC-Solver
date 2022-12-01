/**
 * @file PDivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre 1/Sin D_phi parallel integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivS1Dp.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

#include "QuICC/Transform/Poly/KokkosUtils.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template<>
   PDivS1Dp<CudaIALegendreOperatorTypes>::PDivS1Dp()
   {
   }

   template<>
   PDivS1Dp<CudaIALegendreOperatorTypes>::~PDivS1Dp()
   {
   }

   template <>
   void PDivS1Dp<CudaIALegendreOperatorTypes>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       PDivS1::makeOperator(op, igrid, iweights, i);
       op = static_cast<MHDFloat>(-this->mspSetup->slow(i)) * op;
   }

#ifdef KOKKOS_ENABLE_CUDA
   template <>
   void PDivS1Dp<CudaIALegendreOperatorTypes>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      PDivS1::applyUnitOperator(rOutView, inView, scan, total);

      DataType constant;
      constant.real() = Math::cI.real();
      constant.imag() = Math::cI.imag();

      constantMultiplyMatrix(constant, rOutView);
   }
#endif

   template class PDivS1Dp<CudaIALegendreOperatorTypes>;
}
}
}
}
}
