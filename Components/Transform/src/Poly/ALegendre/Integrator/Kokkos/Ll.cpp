/**
 * @file Ll.cpp
 * @brief Source of the implementation of the associated Legendre P parallel
 * integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/Ll.hpp"
#include "QuICC/Polynomial/ALegendre/Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

void Ll<kokkos_t>::makeOperator(Matrix& op, const Internal::Array& igrid,
   const Internal::Array& iweights, const int i) const
{
   P<kokkos_t>::makeOperator(op, igrid, iweights, i);
   op = op * this->mLl.bottomRows(op.cols()).asDiagonal();
}

void Ll<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   P<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);
}

void Ll<kokkos_t>::initSpecial() const
{
   // Initialise storage for l(l+1) factor
   this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0,
      this->mspSetup->specSize() - 1);
   this->mLl = (this->mLl.array() * (this->mLl.array() + 1.0));
}
} // namespace Integrator
} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC
