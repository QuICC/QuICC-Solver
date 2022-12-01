/**
 * @file PDivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre 1/S D_phi projector
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Projector/PDivS1Dp.hpp"

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

namespace Projector {

   template<>
   PDivS1Dp<CudaIALegendreOperatorTypes>::PDivS1Dp()
   {
   }

   template<>
   PDivS1Dp<CudaIALegendreOperatorTypes>::~PDivS1Dp()
   {
   }

#ifdef KOKKOS_ENABLE_CUDA
   template <>
   void PDivS1Dp<CudaIALegendreOperatorTypes>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      constantMultiplyMatrix<1>(this->mspSetup, scan, inView);
      PDivS1::applyUnitOperator(rOutView, inView, scan, total);
   }
#endif

   template class PDivS1Dp<CudaIALegendreOperatorTypes>;
}
}
}
}
}
