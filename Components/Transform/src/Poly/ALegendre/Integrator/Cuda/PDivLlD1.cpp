/**
 * @file PDivLlD1.cpp
 * @brief Source of the implementation of the associated Legendre l/l(l+1) parallel integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivLlD1.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"

#include "QuICC/Transform/Poly/KokkosUtils.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template<>
   PDivLlD1<CudaIALegendreOperatorTypes>::PDivLlD1()
   {
   }

   template<>
   PDivLlD1<CudaIALegendreOperatorTypes>::~PDivLlD1()
   {
   }

   template <>
   void PDivLlD1<CudaIALegendreOperatorTypes>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       PD1::makeOperator(op, igrid, iweights, i);
       op = op * this->mDivLl.bottomRows(op.cols()).asDiagonal();
   }

#ifdef KOKKOS_ENABLE_CUDA
   template <>
   void PDivLlD1<CudaIALegendreOperatorTypes>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      PD1::applyUnitOperator(rOutView, inView, scan, total);
   }
#endif

   template<>
   void PDivLlD1<CudaIALegendreOperatorTypes>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mDivLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mDivLl = (this->mDivLl.array()*(this->mDivLl.array() + 1.0)).pow(-1);
      this->mDivLl(0) = 0.0;
   }

   template class PDivLlD1<CudaIALegendreOperatorTypes>;
}
}
}
}
}
