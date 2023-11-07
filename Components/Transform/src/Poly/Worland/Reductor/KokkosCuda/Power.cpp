/**
 * @file Power.cpp
 * @brief Source of the implementation of the Worland power spectrum operator
 */


// System includes
//
#include <cassert>


// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/Power.hpp"
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

Power<kokkos_t>::Power() : KokkosIWorlandPower(1)
{
   this->setProfileTag();
}

void Power<kokkos_t>::makeOperator(Matrix& op, Matrix& eop,
   const Internal::Array& igrid, const Internal::Array& iweights,
   const int i) const
{
   int l = this->mspSetup->slow(i);
   int nPoly = this->mspSetup->fastSize(i);

   // Build operator
   op.resize(igrid.size(), nPoly);
   namespace ev = Polynomial::Worland::Evaluator;
   Polynomial::Worland::r_1Wnl<QuICC::Polynomial::Worland::recurrence_t> bwnl;
   bwnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());

   eop.resize(igrid.size(), nPoly);
   Polynomial::Worland::Wnl fwnl(Polynomial::Worland::Wnl::ALPHA_SPHENERGY,
      Polynomial::Worland::Wnl::DBETA_SPHENERGY);
   if (l == 0)
   {
      eop.setZero();
   }
   else
   {
      fwnl.compute<MHDFloat>(eop, nPoly, l - 1, igrid, iweights, ev::Set());
   }
}

void Power<kokkos_t>::applyUnitOperator(const OpMatrixL& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   this->defaultApplyUnitOperator(rOutView, inView, scan, total);
}

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
