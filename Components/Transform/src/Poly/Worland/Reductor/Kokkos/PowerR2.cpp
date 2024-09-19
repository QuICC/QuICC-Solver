/**
 * @file PowerR2.cpp
 * @brief Source of the implementation of the Worland R^2 power spectrum
 * operator
 */


// System includes
//
#include <cassert>


// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/PowerR2.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

PowerR2<kokkos_t>::PowerR2() : KokkosIWorlandPower(1)
{
   this->setProfileTag();
}

void PowerR2<kokkos_t>::makeOperator(Matrix& op, Matrix& eop,
   const Internal::Array& igrid, const Internal::Array& iweights,
   const int i) const
{
   int l = this->mspSetup->slow(i);
   int nPoly = this->mspSetup->fastSize(i);

   // Build operator
   op.resize(igrid.size(), nPoly);
   namespace ev = Polynomial::Worland::Evaluator;
   Polynomial::Worland::Wnl bwnl;
   bwnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());

   Polynomial::Worland::Wnl fwnl(Polynomial::Worland::worland_sphenergy_t::ALPHA,
      Polynomial::Worland::worland_sphenergy_t::DBETA);

   eop.resize(igrid.size(), nPoly);
   fwnl.compute<MHDFloat>(eop, nPoly, l, igrid, iweights, ev::Set());
}

void PowerR2<kokkos_t>::applyUnitOperator(const OpMatrixL& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   this->defaultApplyUnitOperator(rOutView, inView, scan, total);
}

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
