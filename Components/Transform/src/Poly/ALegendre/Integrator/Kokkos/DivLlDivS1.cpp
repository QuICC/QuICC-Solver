/**
 * @file DivLlDivS1.cpp
 * @brief Source of the implementation of the associated Legendre 1/l(l+1) 1/Sin
 * parallel integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/DivLlDivS1.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

void DivLlDivS1<kokkos_t>::makeOperator(Matrix& op,
   const Internal::Array& igrid, const Internal::Array& iweights,
   const int i) const
{
   DivS1<kokkos_t>::makeOperator(op, igrid, iweights, i);
   op = op * this->mDivLl.bottomRows(op.cols()).asDiagonal();
}

void DivLlDivS1<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   DivS1<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);
}

void DivLlDivS1<kokkos_t>::initSpecial() const
{
   // Initialise storage for l(l+1) factor
   this->mDivLl = Array::LinSpaced(this->mspSetup->specSize(), 0,
      this->mspSetup->specSize() - 1);
   this->mDivLl = (this->mDivLl.array() * (this->mDivLl.array() + 1.0)).pow(-1);
   this->mDivLl(0) = 0.0;
}

} // namespace Integrator
} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC
