/**
 * @file I3_Zero.cpp
 * @brief Source of the implementation of the associated Worland I3 parallel
 * integrator
 */


// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/I3_Zero.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/SparseSM/Worland/I3.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

I3_Zero<kokkos_t>::I3_Zero()
{
   this->setProfileTag();
}

void I3_Zero<kokkos_t>::makeOperator(Matrix& op, const Internal::Array& igrid,
   const Internal::Array& iweights, const int i) const
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
      // Internal computation uses dealiased modes
      const int extraN = 5 * (!this->mcTruncQI); // I3 has 5 superdiagonals
      int nN = nPoly + extraN;
      this->checkGridSize(nN, l, igrid.size());

      Internal::Matrix tOp(igrid.size(), nN);

      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;
      wnl.compute<Internal::MHDFloat>(tOp, nN, l, igrid, iweights, ev::Set());

      auto a = wnl.alpha(l);
      auto b = wnl.dBeta();
      ::QuICC::SparseSM::Worland::I3 spasm(nN, nN, a, b, l,
         1 * this->mcTruncQI);
      tOp = (spasm.mat() * tOp.transpose()).transpose();
      op = tOp.cast<MHDFloat>().leftCols(nPoly);
   }
}

void I3_Zero<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   applyBlockOperator<5>(this->mspSetup, this->vmOps, rOutView, inView, scan,
      total);
}

} // namespace Integrator
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
