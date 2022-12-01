/**
 * @file PDivLlDivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre 1/l(l+1) 1/Sin D_phi parallel integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivLlDivS1Dp.hpp"

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
   PDivLlDivS1Dp<CudaIALegendreOperatorTypes>::PDivLlDivS1Dp()
   {
   }

   template<>
   PDivLlDivS1Dp<CudaIALegendreOperatorTypes>::~PDivLlDivS1Dp()
   {
   }

   template <>
   void PDivLlDivS1Dp<CudaIALegendreOperatorTypes>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       PDivS1::makeOperator(op, igrid, iweights, i);
       op = op * this->mDivLl.bottomRows(op.cols()).asDiagonal();
   }

#ifdef KOKKOS_ENABLE_CUDA
   template <>
   void PDivLlDivS1Dp<CudaIALegendreOperatorTypes>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      PDivS1::applyUnitOperator(rOutView, inView, scan, total);
      constantMultiplyMatrix(this->mspSetup, scan, rOutView);
   }
#endif

   template<>
   void PDivLlDivS1Dp<CudaIALegendreOperatorTypes>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mDivLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mDivLl = (this->mDivLl.array()*(this->mDivLl.array() + 1.0)).pow(-1);
      this->mDivLl(0) = 0.0;
   }

   template class PDivLlDivS1Dp<CudaIALegendreOperatorTypes>;
}
}
}
}
}
