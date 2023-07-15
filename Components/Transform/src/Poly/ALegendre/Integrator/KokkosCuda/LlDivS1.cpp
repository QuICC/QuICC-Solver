/**
 * @file LlDivS1.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1) 1/Sin parallel integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/LlDivS1.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

#include "QuICC/Transform/Poly/KokkosUtils.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   void LlDivS1<kokkos_t>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       DivS1<kokkos_t>::makeOperator(op, igrid, iweights, i);
       op = op * this->mLl.bottomRows(op.cols()).asDiagonal();
   }

   void LlDivS1<kokkos_t>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      DivS1<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);
   }

   void LlDivS1<kokkos_t>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl = (this->mLl.array()*(this->mLl.array() + 1.0));
   }

}
}
}
}
}
