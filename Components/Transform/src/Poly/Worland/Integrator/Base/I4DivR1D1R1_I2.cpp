/**
 * @file I4DivR1D1R1_I2.cpp
 * @brief Source of the implementation of the Worland 1/R1 D R1 integrator but 0 mode is I2 P integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/I4DivR1D1R1_I2.hpp"
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

   I4DivR1D1R1_I2<base_t>::I4DivR1D1R1_I2()
   {
   }

   void I4DivR1D1R1_I2<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      if(l == 0)
      {
         l = 1;

         // Internal computation uses dealiased modes
         const int extraN = 3*(!this->mcTruncQI); // I2 has 3 superdiagonals
         int nN = nPoly + extraN;
         this->checkGridSize(nN, l, igrid.size());

         Internal::Matrix tOp(igrid.size(), nN);

         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;

         wnl.compute<Internal::MHDFloat>(tOp, nN, l, igrid, iweights, ev::Set());

         auto a = wnl.alpha(l);
         auto b = wnl.dBeta();
         ::QuICC::SparseSM::Worland::I2 spasm(nN, nN, a, b, l, 1*this->mcTruncQI);
         tOp = (spasm.mat()*tOp.transpose()).transpose();
         op = tOp.cast<MHDFloat>().leftCols(nPoly);
      } else
      {
         // Internal computation uses dealiased modes
         const int extraN = 6*(!this->mcTruncQI); // I4 has 6 superdiagonals
         int nN = nPoly + extraN;
         this->checkGridSize(nN, l, igrid.size());

         Internal::Matrix tOp(igrid.size(), nN);

         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;

#if defined QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
         // **************************************************
         // Formulation without explicit grid:
         // Operates on polynomials with l = l-1
         int l_in = std::abs(l-1);
         int n_in = nN + 1;
         this->checkGridSize(n_in, l_in, igrid.size());

         Internal::Matrix opA(igrid.size(), n_in);
         wnl.compute<Internal::MHDFloat>(opA, n_in, l_in, igrid, iweights, ev::Set());

         Internal::Matrix opB(igrid.size(), n_in);
         Polynomial::Worland::r_1drWnl r_1drWnl;
         r_1drWnl.compute<Internal::MHDFloat>(opB, n_in, l_in, igrid, Internal::Array(), ev::Set());

         Internal::Matrix opC(igrid.size(), nN);
         Polynomial::Worland::Wnl wnlB;
         wnlB.compute<Internal::MHDFloat>(opC, nN, l, igrid, iweights, ev::Set());

         tOp = (opC.transpose()*opB*opA.transpose()).transpose();
#else

         // **************************************************
         // Alternative formulation of operators:
         // This version uses explicit radial factor to work on l polynomials

         Internal::Matrix opA(igrid.size(), nN);
         wnl.compute<Internal::MHDFloat>(opA, nN, l, igrid, iweights, ev::Set());

         Internal::Matrix opB(igrid.size(), nN);
         Polynomial::Worland::dWnl dWnl;
         dWnl.compute<Internal::MHDFloat>(opB, nN, l, igrid, Internal::Array(), ev::Set());

         tOp = (opA.transpose()*igrid.array().pow(-1).matrix().asDiagonal()*opB*opA.transpose()*igrid.asDiagonal()).transpose();
#endif

         // Multiply by Quasi-inverse
         auto a = wnl.alpha(l);
         auto b = wnl.dBeta();
         ::QuICC::SparseSM::Worland::I4 spasm(nN, nN, a, b, l, 2*this->mcTruncQI);
         tOp = (spasm.mat()*tOp.transpose()).transpose();
         op = tOp.cast<MHDFloat>().leftCols(nPoly);

         assert(op.rows() == igrid.size());
         assert(op.cols() == nPoly);
      }
   }

   void I4DivR1D1R1_I2<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}

#undef QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
