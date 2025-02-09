/**
 * @file DivLlDivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre 1/l(l+1) 1/Sin D_phi parallel integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/DivLlDivS1Dp.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"
#include "QuICC/Transform/Poly/KokkosUtils.hpp"
#include "QuICC/Transform/Poly/ALegendre/CudaIALegendreOperatorGemmUtils.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   void DivLlDivS1Dp<kokkos_t>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       DivS1<kokkos_t>::makeOperator(op, igrid, iweights, i);
       op = op * this->mDivLl.bottomRows(op.cols()).asDiagonal();
   }

   void DivLlDivS1Dp<kokkos_t>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      DivS1<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);
      constantMultiplyMatrix(this->mspSetup, scan, rOutView);
   }

   void DivLlDivS1Dp<kokkos_t>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mDivLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mDivLl = (this->mDivLl.array()*(this->mDivLl.array() + 1.0)).pow(-1);
      this->mDivLl(0) = 0.0;
   }

}
}
}
}
}
