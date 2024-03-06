/**
 * @file DivLlD1.cpp
 * @brief Source of the implementation of the associated Legendre l/l(l+1)
 * parallel integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/DivLlD1.hpp"
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

void DivLlD1<kokkos_t>::makeOperator(Matrix& op, const Internal::Array& igrid,
   const Internal::Array& iweights, const int i) const
{
   D1<kokkos_t>::makeOperator(op, igrid, iweights, i);
   op = op * this->mDivLl.bottomRows(op.cols()).asDiagonal();
}

void DivLlD1<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   D1<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);
}

void DivLlD1<kokkos_t>::initSpecial() const
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
