/**
 * @file PLlDivS1.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1) 1/Sin parallel integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/PLlDivS1.hpp"

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
   PLlDivS1<CudaIALegendreOperatorTypes>::PLlDivS1()
   {
   }

   template<>
   PLlDivS1<CudaIALegendreOperatorTypes>::~PLlDivS1()
   {
   }

   template <>
   void PLlDivS1<CudaIALegendreOperatorTypes>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       PDivS1::makeOperator(op, igrid, iweights, i);
       op = op * this->mLl.bottomRows(op.cols()).asDiagonal();
   }

#ifdef KOKKOS_ENABLE_CUDA
   template <>
   void PLlDivS1<CudaIALegendreOperatorTypes>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      PDivS1::applyUnitOperator(rOutView, inView, scan, total);
   }
#endif

   template<>
   void PLlDivS1<CudaIALegendreOperatorTypes>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl = (this->mLl.array()*(this->mLl.array() + 1.0));
   }

   template class PLlDivS1<CudaIALegendreOperatorTypes>;
}
}
}
}
}
