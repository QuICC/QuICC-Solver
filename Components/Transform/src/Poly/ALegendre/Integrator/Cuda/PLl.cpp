/**
 * @file PLl.cpp
 * @brief Source of the implementation of the associated Legendre P parallel integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/PLl.hpp"

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
   PLl<CudaIALegendreOperatorTypes>::PLl()
   {
   }

   template<>
   PLl<CudaIALegendreOperatorTypes>::~PLl()
   {
   }

   template <>
   void PLl<CudaIALegendreOperatorTypes>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       PP::makeOperator(op, igrid, iweights, i);
       op = op * this->mLl.bottomRows(op.cols()).asDiagonal();
   }

   template<>
   void PLl<CudaIALegendreOperatorTypes>::applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const
   {
      #if defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
      rOut = this->mOps.at(i).transpose() * in;
      #endif //defined QUICC_ALEGENDRE_INTGIMPL_OTF
   }

#ifdef KOKKOS_ENABLE_CUDA
   template <>
   void PLl<CudaIALegendreOperatorTypes>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      PP::applyUnitOperator(rOutView, inView, scan, total);
   }
#endif

   template<>
   void PLl<CudaIALegendreOperatorTypes>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl = (this->mLl.array()*(this->mLl.array() + 1.0));
   }

   template class PLl<CudaIALegendreOperatorTypes>;
}
}
}
}
}
