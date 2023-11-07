/**
 * @file D1R1.cpp
 * @brief Source of the implementation of the Worland D R projector
 */

// System includes
//
#include <cassert>


// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/D1R1.hpp"
#include "QuICC/Polynomial/Worland/drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

D1R1<kokkos_t>::D1R1()
{
   this->setProfileTag();
}

void D1R1<kokkos_t>::makeOperator(Matrix& op, const Internal::Array& igrid,
   const Internal::Array& iweights, const int i) const
{
   int l = this->mspSetup->slow(i);

   // Build operator
   int nPoly = this->mspSetup->fastSize(i);
   op.resize(igrid.size(), nPoly);
   namespace ev = Polynomial::Worland::Evaluator;
   Polynomial::Worland::drWnl wnl;
   wnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());
}

void D1R1<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   applyBlockOperator<4>(this->mspSetup, this->vmOps, rOutView, inView, scan,
      total);
}

} // namespace Projector
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
