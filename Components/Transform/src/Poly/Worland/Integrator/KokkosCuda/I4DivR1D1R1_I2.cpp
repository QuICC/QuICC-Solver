/**
 * @file I4DivR1D1R1_I2.cpp
 * @brief Source of the implementation of the associated Worland I2 1/R D R
 * parallel integrator but 0 mode zeroed.
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/I4DivR1D1R1_I2.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/dWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"

// It is not clear yet which implementation is more accurate
#define QUICC_AVOID_EXPLICIT_RADIAL_FACTOR

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

void I4DivR1D1R1_I2<kokkos_t>::makeOperator(Matrix& op,
   const Internal::Array& igrid, const Internal::Array& iweights,
   const int i) const
{
   int l = this->mspSetup->slow(i);

   // Build operator
   int nPoly = this->mspSetup->fastSize(i);
   if (l == 0)
   {
      l = 1;

      // Internal computation uses dealiased modes
      const int extraN = 3 * (!this->mcTruncQI); // I2 has 3 superdiagonals
      int nN = nPoly + extraN;
      this->checkGridSize(nN, l, igrid.size());

      Internal::Matrix tOp(igrid.size(), nN);

      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;

      wnl.compute<Internal::MHDFloat>(tOp, nN, l, igrid, iweights, ev::Set());

      auto a = wnl.alpha(l);
      auto b = wnl.dBeta();
      ::QuICC::SparseSM::Worland::I2 spasm(nN, nN, a, b, l,
         1 * this->mcTruncQI);
      tOp = (spasm.mat() * tOp.transpose()).transpose();
      op = tOp.cast<MHDFloat>().leftCols(nPoly);
   }
   else
   {
      // Internal computation uses dealiased modes
      const int extraN = 6 * (!this->mcTruncQI); // I4 has 6 superdiagonals
      int nN = nPoly + extraN;
      this->checkGridSize(nN, l, igrid.size());

      Internal::Matrix tOp(igrid.size(), nN);

      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;

#ifdef QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
      using poly_t = QuICC::Polynomial::Worland::implicit_t;
#else
      using poly_t = QuICC::Polynomial::Worland::explicit_t;
#endif
      Polynomial::Worland::r_1drWnl<poly_t> r_1drWnl;
      r_1drWnl.compute<Internal::MHDFloat>(tOp, nN, l, igrid, iweights,
         ev::Set());


      // Multiply by Quasi-inverse
      auto a = wnl.alpha(l);
      auto b = wnl.dBeta();
      ::QuICC::SparseSM::Worland::I4 spasm(nN, nN, a, b, l,
         2 * this->mcTruncQI);
      tOp = (spasm.mat() * tOp.transpose()).transpose();
      op = tOp.cast<MHDFloat>().leftCols(nPoly);

      assert(op.rows() == igrid.size());
      assert(op.cols() == nPoly);
   }
}

void I4DivR1D1R1_I2<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
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
