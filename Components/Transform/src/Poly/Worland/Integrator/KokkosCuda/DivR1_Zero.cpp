/**
 * @file DivR1_Zero.cpp
 * @brief Source of the implementation of the Worland 1/R integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/DivR1_Zero.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

// It is not clear yet which implementation is more accurate
#define QUICC_AVOID_EXPLICIT_RADIAL_FACTOR


namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

DivR1_Zero<kokkos_t>::DivR1_Zero()
{
   this->setProfileTag();
}

void DivR1_Zero<kokkos_t>::makeOperator(Matrix& op,
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
#if defined QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
      using poly_t = QuICC::Polynomial::Worland::implicit_t;
#else
      using poly_t = QuICC::Polynomial::Worland::explicit_t;
#endif
      Polynomial::Worland::r_1Wnl<poly_t> wnl;
      wnl.compute<MHDFloat>(op, nPoly, l, igrid, iweights, ev::Set());
   }
}

void DivR1_Zero<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   applyBlockOperator<3>(this->mspSetup, this->vmOps, rOutView, inView, scan,
      total);
}

} // namespace Integrator
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#undef QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
