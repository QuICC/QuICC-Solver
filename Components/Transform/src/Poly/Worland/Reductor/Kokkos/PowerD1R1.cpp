/**
 * @file PowerD1R1.cpp
 * @brief Source of the implementation of the Worland D R power spectrum
 * operator
 */


// System includes
//
#include <cassert>


// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/PowerD1R1.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"



namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

PowerD1R1<kokkos_t>::PowerD1R1() : KokkosIWorlandPower(0)
{
   this->setProfileTag();
}

void PowerD1R1<kokkos_t>::makeOperator(Matrix& op, Matrix& eop,
   const Internal::Array& igrid, const Internal::Array& iweights,
   const int i) const
{
   int l = this->mspSetup->slow(i);
   int nPoly = this->mspSetup->fastSize(i);

   // Build operator
   op.resize(igrid.size(), nPoly);
   namespace ev = Polynomial::Worland::Evaluator;
   Polynomial::Worland::r_1drWnl<QuICC::Polynomial::Worland::recurrence_t> bwnl;
   bwnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());

   eop.resize(igrid.size(), nPoly);
   Polynomial::Worland::Wnl fwnl(Polynomial::Worland::worland_sphenergy_t::ALPHA,
      Polynomial::Worland::worland_sphenergy_t::DBETA);
   if (l == 0)
   {
      eop.setZero();
   }
   else
   {
      fwnl.compute<MHDFloat>(eop, nPoly, std::abs(l - 1), igrid, iweights,
         ev::Set());
   }
}

void PowerD1R1<kokkos_t>::applyUnitOperator(const OpMatrixL& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   this->defaultApplyUnitOperator(rOutView, inView, scan, total);
}

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
