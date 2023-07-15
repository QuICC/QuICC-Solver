/**
 * @file DivLl.cpp
 * @brief Source of the implementation of the associated Legendre l/l(l+1)  parallel integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/DivLl.hpp"
#include "QuICC/Polynomial/ALegendre/Plm.hpp"

#include "QuICC/Transform/Poly/KokkosUtils.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
/* #include <type_traits> */

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   void DivLl<kokkos_t>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
       P<kokkos_t>::makeOperator(op, igrid, iweights, i);
       op = op * this->mDivLl.bottomRows(op.cols()).asDiagonal();
   }

   void DivLl<kokkos_t>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      P<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);
   }

   void DivLl<kokkos_t>::initSpecial() const
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
