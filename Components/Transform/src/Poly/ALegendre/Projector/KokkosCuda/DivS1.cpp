/**
 * @file DivS1.cpp
 * @brief Source of the implementation of the associated Legendre 1/Sin
 * projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/DivS1.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

void DivS1<kokkos_t>::makeOperator(Matrix& op, const Internal::Array& igrid,
   const Internal::Array& iweights, const int i) const
{
   int m = this->mspSetup->slow(i);
   int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i) - 1, i) - m + 1;

   // Build operator
   op.resize(igrid.size(), nPoly);
   namespace ev = Polynomial::ALegendre::Evaluator;
   Polynomial::ALegendre::sin_1Plm plm;
   plm.compute<MHDFloat>(op, nPoly, m, igrid, Internal::Array(), ev::Set());
}

void DivS1<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   applyBlockOperator<1>(this->mspSetup, this->vmOps, rOutView, inView, scan,
      total);
}

} // namespace Projector
} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC
