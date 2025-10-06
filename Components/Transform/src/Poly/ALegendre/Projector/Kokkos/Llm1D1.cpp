/**
 * @file Llm1D1.cpp
 * @brief Source of the parallel implementation of the associated Legendre
 * [l(l+1)-1] D projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/Llm1D1.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

void Llm1D1<kokkos_t>::makeOperator(Matrix& op, const Internal::Array& igrid,
   const Internal::Array& iweights, const int i) const
{
   D1<kokkos_t>::makeOperator(op, igrid, iweights, i);
   op = op * this->mLlm1.bottomRows(op.cols()).asDiagonal();
}

void Llm1D1<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   D1<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);
}

void Llm1D1<kokkos_t>::initSpecial() const
{
   // Initialise storage for [l(l+1)-1] factor
   this->mLlm1 = Array::LinSpaced(this->mspSetup->specSize(), 0,
      this->mspSetup->specSize() - 1);
   this->mLlm1 = this->mLlm1.array() * (this->mLlm1.array() + 1.0)-1.0;
}

} // namespace Projector
} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC
