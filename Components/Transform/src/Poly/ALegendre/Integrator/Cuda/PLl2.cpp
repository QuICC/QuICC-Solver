/**
 * @file PLl2.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1)^2 P parallel integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/PLl2.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/Plm.hpp"

#include "QuICC/Transform/Poly/KokkosUtils.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template<>
   PLl2<CudaIALegendreOperatorTypes>::PLl2()
   {
   }

   template<>
   PLl2<CudaIALegendreOperatorTypes>::~PLl2()
   {
   }

   template <>
   void PLl2<CudaIALegendreOperatorTypes>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       PP::makeOperator(op, igrid, iweights, i);
       op = op * this->mLl2.bottomRows(op.cols()).asDiagonal();
   }

   template<>
   void PLl2<CudaIALegendreOperatorTypes>::applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const
   {
      #if defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
      rOut = this->mOps.at(i).transpose() * in;
      #endif //defined QUICC_ALEGENDRE_INTGIMPL_OTF
   }

#ifdef KOKKOS_ENABLE_CUDA
   template <>
   void PLl2<CudaIALegendreOperatorTypes>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      PP::applyUnitOperator(rOutView, inView, scan, total);
   }
#endif

   template<>
   void PLl2<CudaIALegendreOperatorTypes>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLl2 = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl2 = (this->mLl2.array()*(this->mLl2.array() + 1.0)).pow(2);
   }

   template class PLl2<CudaIALegendreOperatorTypes>;
}
}
}
}
}
