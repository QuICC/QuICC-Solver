/**
 * @file PDivLlDivS1.cpp
 * @brief Source of the implementation of the associated Legendre 1/l(l+1) 1/Sin parallel integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivLlDivS1.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

#include "QuICC/Transform/Poly/KokkosUtils.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template<>
   PDivLlDivS1<CudaIALegendreOperatorTypes>::PDivLlDivS1()
   {
   }

   template<>
   PDivLlDivS1<CudaIALegendreOperatorTypes>::~PDivLlDivS1()
   {
   }

   template <>
   void PDivLlDivS1<CudaIALegendreOperatorTypes>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       PDivS1::makeOperator(op, igrid, iweights, i);
       op = op * this->mDivLl.bottomRows(op.cols()).asDiagonal();
   }

#ifdef KOKKOS_ENABLE_CUDA
   template <>
   void PDivLlDivS1<CudaIALegendreOperatorTypes>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      PDivS1::applyUnitOperator(rOutView, inView, scan, total);
   }
#endif

   template<>
   void PDivLlDivS1<CudaIALegendreOperatorTypes>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mDivLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mDivLl = (this->mDivLl.array()*(this->mDivLl.array() + 1.0)).pow(-1);
      this->mDivLl(0) = 0.0;
   }

   template class PDivLlDivS1<CudaIALegendreOperatorTypes>;
}
}
}
}
}
