/**
 * @file DivR1_Zero.cpp
 * @brief Source of the implementation of the Worland 1/R integrator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/DivR1_Zero.hpp"
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

   DivR1_Zero<base_t>::DivR1_Zero()
   {
      this->setProfileTag();
   }

   void DivR1_Zero<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      if(l == 0)
      {
         op.resize(igrid.size(), nPoly);
         op.setZero();
      } else
      {
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;

         // Internal computation uses dealiased modes
         int nN = nPoly + 0;
         this->checkGridSize(nN, l, igrid.size());

         Internal::Matrix tOp(igrid.size(), nN);
#if defined QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
         // **************************************************
         // Formulation without explicit grid:
         // Operates on polynomials with l = l-1
         int l_in = std::abs(l-1);
         int l_out = std::abs(l-0);
         int n_in = nN + 1;
         this->checkGridSize(n_in, l_in, igrid.size());

         Internal::Matrix opA(igrid.size(), n_in);
         wnl.compute<Internal::MHDFloat>(opA, n_in, l_in, igrid, iweights, ev::Set());

         Internal::Matrix opB(igrid.size(), n_in);
         Polynomial::Worland::r_1Wnl r_1Wnl;
         r_1Wnl.compute<Internal::MHDFloat>(opB, n_in, l_in, igrid, Internal::Array(), ev::Set());

         Internal::Matrix opC(igrid.size(), nN);
         Polynomial::Worland::Wnl wnlB;
         wnlB.compute<Internal::MHDFloat>(opC, nN, l_out, igrid, iweights, ev::Set());

         tOp = (opC.transpose()*opB*opA.transpose()).transpose();
#else

         // **************************************************
         // Alternative formulation of operators:
         // This version uses explicit radial factors to work on l polynomials

         Internal::Matrix opA(igrid.size(), nN);
         wnl.compute<Internal::MHDFloat>(opA, nN, l, igrid, iweights, ev::Set());

         tOp = (opA.transpose()*igrid.array().pow(-1).matrix().asDiagonal()).transpose();
#endif
         op = tOp.cast<MHDFloat>().leftCols(nPoly);

         assert(op.rows() == igrid.size());
         assert(op.cols() == nPoly);
      }
   }

   void DivR1_Zero<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}

#undef QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
