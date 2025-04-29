/**
 * @file P.cpp
 * @brief Source of the implementation of the associated Legendre P parallel
 * integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/P.hpp"
#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/InnerProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

void P<kokkos_t>::makeOperator(Matrix& op, const Internal::Array& igrid,
   const Internal::Array& iweights, const int i) const
{
   int m = this->mspSetup->slow(i);
   int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i) - 1, i) - m + 1;

   // Build operator
   op.resize(igrid.size(), nPoly);
   namespace ev = Polynomial::ALegendre::Evaluator;
   Polynomial::ALegendre::Plm plm;
   plm.compute<MHDFloat>(op, nPoly, m, igrid, iweights, ev::Set());
}

void P<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   applyBlockOperator(this->mspSetup, this->vmOps, rOutView, inView, scan,
      total);
}

} // namespace Integrator
} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC
