/**
 * @file Ll2.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1)^2 P
 * parallel integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/Ll2.hpp"
#include "QuICC/Polynomial/ALegendre/Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

void Ll2<kokkos_t>::makeOperator(Matrix& op, const Internal::Array& igrid,
   const Internal::Array& iweights, const int i) const
{
   P<kokkos_t>::makeOperator(op, igrid, iweights, i);
   op = op * this->mLl2.bottomRows(op.cols()).asDiagonal();
}

void Ll2<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   P<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);
}

void Ll2<kokkos_t>::initSpecial() const
{
   // Initialise storage for l(l+1) factor
   this->mLl2 = Array::LinSpaced(this->mspSetup->specSize(), 0,
      this->mspSetup->specSize() - 1);
   this->mLl2 = (this->mLl2.array() * (this->mLl2.array() + 1.0)).pow(2);
}

} // namespace Integrator
} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC
