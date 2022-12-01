/**
 * @file PDivLl.cpp
 * @brief Source of the implementation of the associated Legendre l/l(l+1)  parallel integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivLl.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/Plm.hpp"

#include "QuICC/Transform/Poly/KokkosUtils.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
/* #include <type_traits> */

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template<>
   PDivLl<CudaIALegendreOperatorTypes>::PDivLl()
   {
   }

   template<>
   PDivLl<CudaIALegendreOperatorTypes>::~PDivLl()
   {
   }

   template <>
   void PDivLl<CudaIALegendreOperatorTypes>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       PP::makeOperator(op, igrid, iweights, i);
       op = op * this->mDivLl.bottomRows(op.cols()).asDiagonal();
   }

   template<>
   void PDivLl<CudaIALegendreOperatorTypes>::applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const
   {
      #if defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
      rOut = this->mOps.at(i).transpose() * in;
      #endif //defined QUICC_ALEGENDRE_INTGIMPL_OTF
   }

#ifdef KOKKOS_ENABLE_CUDA
   template <>
   void PDivLl<CudaIALegendreOperatorTypes>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      PP::applyUnitOperator(rOutView, inView, scan, total);
   }
#endif

   template<>
   void PDivLl<CudaIALegendreOperatorTypes>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mDivLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mDivLl = (this->mDivLl.array()*(this->mDivLl.array() + 1.0)).pow(-1);
      this->mDivLl(0) = 0.0;
   }

   template class PDivLl<CudaIALegendreOperatorTypes>;
}
}
}
}
}
