/**
 * @file RadialPowerDivR1D1R1.cpp
 * @brief Source of the implementation of the Worland 1/R D R radial power
 * spectrum operator
 */


// System includes
//
#include <cassert>


// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/RadialPowerDivR1D1R1.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

RadialPowerDivR1D1R1<kokkos_t>::RadialPowerDivR1D1R1() :
    KokkosIWorlandRadialPower()
{
   this->setProfileTag();
}

void RadialPowerDivR1D1R1<kokkos_t>::makeOperator(Matrix& op,
   const Internal::Array& igrid, const Internal::Array& iweights,
   const int i) const
{
   int l = this->mspSetup->slow(i);

   // Build operator
   int nPoly = this->mspSetup->fastSize(i);
   op.resize(igrid.size(), nPoly);
   if (l == 0)
   {
      op.setZero();
   }
   else
   {
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::r_1drWnl<QuICC::Polynomial::Worland::recurrence_t>
         wnl;
      wnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());
   }
}

void RadialPowerDivR1D1R1<kokkos_t>::applyUnitOperator(
   const OpMatrixL& rOutView, const OpMatrixLZ& inView, const OpVectorI& scan,
   const int total) const
{
   applyKokkosBlockOperator<4>(this->mspSetup, this->vmOps, rOutView, inView,
      scan, Abs2Complex());
}

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
