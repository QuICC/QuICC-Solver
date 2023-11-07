/**
 * @file D1CylLaplh.cpp
 * @brief Source of the implementation of the Worland D of cylindrical
 * horizontal laplacian projector
 */


// System includes
//
#include <cassert>


// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/D1CylLaplh.hpp"
#include "QuICC/Polynomial/Worland/dclaplhWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

void D1CylLaplh<kokkos_t>::makeOperator(Matrix& op,
   const Internal::Array& igrid, const Internal::Array& iweights,
   const int i) const
{
   int l = this->mspSetup->slow(i);

   // Build operator
   int nPoly = this->mspSetup->fastSize(i);
   op.resize(igrid.size(), nPoly);
   namespace ev = Polynomial::Worland::Evaluator;
   Polynomial::Worland::dclaplhWnl wnl;
   wnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());
}

void D1CylLaplh<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
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
