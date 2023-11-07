/**
 * @file RadialPower.cpp
 * @brief Source of the implementation of the Worland power spectrum operator
 */


// System includes
//
#include <cassert>


// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/RadialPower.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

RadialPower<kokkos_t>::RadialPower() : KokkosIWorlandRadialPower()
{
   this->setProfileTag();
}

void RadialPower<kokkos_t>::makeOperator(Matrix& op,
   const Internal::Array& igrid, const Internal::Array& iweights,
   const int i) const
{
   int l = this->mspSetup->slow(i);

   // Build operator
   int nPoly = this->mspSetup->fastSize(i);
   op.resize(igrid.size(), nPoly);
   namespace ev = Polynomial::Worland::Evaluator;
   Polynomial::Worland::Wnl wnl;
   wnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());
}

void RadialPower<kokkos_t>::applyUnitOperator(const OpMatrixL& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   // Changing col size or row size? If change then it behaves as a worland
   // projector
   auto change = this->mspSetup->mult(this->mspSetup->slowSize() - 1) -
                 this->mspSetup->mult(0);
   if (change)
      applyBlockOperator<4>(this->mspSetup, this->vmOps, rOutView, inView, scan,
         total, Abs2Complex());
   else
      applyBlockOperator<1>(this->mspSetup, this->vmOps, rOutView, inView, scan,
         total, Abs2Complex());
}

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
