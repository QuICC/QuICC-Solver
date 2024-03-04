/**
 * @file LlD1.cpp
 * @brief Source of the parallel implementation of the associated Legendre
 * l(l+1) D projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/LlD1.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

void LlD1<kokkos_t>::makeOperator(Matrix& op, const Internal::Array& igrid,
   const Internal::Array& iweights, const int i) const
{
   D1<kokkos_t>::makeOperator(op, igrid, iweights, i);
   op = op * this->mLl.bottomRows(op.cols()).asDiagonal();
}

void LlD1<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   D1<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);
}

void LlD1<kokkos_t>::initSpecial() const
{
   // Initialise storage for l(l+1) factor
   this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0,
      this->mspSetup->specSize() - 1);
   this->mLl = this->mLl.array() * (this->mLl.array() + 1.0);
}

} // namespace Projector
} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC
