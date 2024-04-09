/**
 * @file P_Zero_Zero.cpp
 * @brief Source of the implementation of the associated Worland P_Zero parallel
 * integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/P_Zero.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"

#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

P_Zero<kokkos_t>::P_Zero() : KokkosIWorlandIntegrator()
{
   this->setProfileTag();
}


void P_Zero<kokkos_t>::makeOperator(Matrix& op, const Internal::Array& igrid,
   const Internal::Array& iweights, const int i) const
{
   int l = this->mspSetup->slow(i);

   // Build operator
   int nPoly = this->mspSetup->fastSize(i);
   if (l == 0)
   {
      op.resize(igrid.size(), nPoly);
      op.setZero();
   }
   else
   {
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;

      // Internal computation uses dealiased modes
      int nN = nPoly;
      this->checkGridSize(nN, l, igrid.size());

      Internal::Matrix tOp(igrid.size(), nN);

      wnl.compute<Internal::MHDFloat>(tOp, nN, l, igrid, iweights, ev::Set());

      op = tOp.cast<MHDFloat>().leftCols(nPoly);
   }
}

void P_Zero<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   applyKokkosBlockOperator<3>(this->mspSetup, this->vmOps, rOutView, inView,
      scan);
}
} // namespace Integrator
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
