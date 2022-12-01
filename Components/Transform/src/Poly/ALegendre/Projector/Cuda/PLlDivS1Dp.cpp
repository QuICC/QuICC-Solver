/**
 * @file PLlDivS1Dp.cpp
 * @brief Source of the parallel implementation of the associated Legendre l(l+1)/Sin D_phi projector
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Projector/PLlDivS1Dp.hpp"

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
   PLlDivS1Dp<CudaIALegendreOperatorTypes>::PLlDivS1Dp()
   {
   }

   template<>
   PLlDivS1Dp<CudaIALegendreOperatorTypes>::~PLlDivS1Dp()
   {
   }

   template <>
   void PLlDivS1Dp<CudaIALegendreOperatorTypes>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       PDivS1::makeOperator(op, igrid, iweights, i);
       op = op * (static_cast<MHDFloat>(this->mspSetup->slow(i)) * this->mLl.bottomRows(op.cols())).asDiagonal();
   }

#ifdef KOKKOS_ENABLE_CUDA
   template <>
   void PLlDivS1Dp<CudaIALegendreOperatorTypes>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      PDivS1::applyUnitOperator(rOutView, inView, scan, total);

      DataType constant;
      constant.real() = Math::cI.real();
      constant.imag() = Math::cI.imag();

      constantMultiplyMatrix(constant, rOutView);
   }
#endif

   template<>
   void PLlDivS1Dp<CudaIALegendreOperatorTypes>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl = this->mLl.array()*(this->mLl.array() + 1.0);
   }

   template class PLlDivS1Dp<CudaIALegendreOperatorTypes>;
}
}
}
}
}
