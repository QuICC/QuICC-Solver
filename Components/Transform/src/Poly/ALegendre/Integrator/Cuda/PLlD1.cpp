/**
 * @file PLlD1.cpp
 * @brief Source of the implementation of the associated Legendre P parallel integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/PLlD1.hpp"

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
   PLlD1<CudaIALegendreOperatorTypes>::PLlD1()
   {
   }

   template<>
   PLlD1<CudaIALegendreOperatorTypes>::~PLlD1()
   {
   }

   template <>
   void PLlD1<CudaIALegendreOperatorTypes>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       PD1::makeOperator(op, igrid, iweights, i);
       op = op * this->mLl.bottomRows(op.cols()).asDiagonal();
   }

#ifdef KOKKOS_ENABLE_CUDA
   template <>
   void PLlD1<CudaIALegendreOperatorTypes>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      PD1::applyUnitOperator(rOutView, inView, scan, total);
   }
#endif

   template<>
   void PLlD1<CudaIALegendreOperatorTypes>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl = (this->mLl.array()*(this->mLl.array() + 1.0));
   }

   template class PLlD1<CudaIALegendreOperatorTypes>;
}
}
}
}
}
