/**
 * @file DivR1D1R1.cpp
 * @brief Source of the implementation of the Worland 1/R D R integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Integrator/DivR1D1R1.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/dWnl.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

// It is not clear yet which implementation is more accurate
#undef QUICC_AVOID_EXPLICIT_RADIAL_FACTOR

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   DivR1D1R1::DivR1D1R1()
   {
      this->setProfileTag();
   }

   DivR1D1R1::~DivR1D1R1()
   {
   }

   void DivR1D1R1::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);

      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;

      // Internal computation uses dealiased modes
      int nN = nPoly + 0;
      this->checkGridSize(nN, l, igrid.size());

      Internal::Matrix tOp(igrid.size(), nN);

#ifdef QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
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
      wnl.compute<Internal::MHDFloat>(opA, nN, l, igrid, iweights.array()*igrid.array(), ev::Set());

      Internal::Matrix opB(igrid.size(), nN);
      Polynomial::Worland::dWnl dWnl;
      dWnl.compute<Internal::MHDFloat>(opB, nN, l, igrid, Internal::Array(), ev::Set());

      Internal::Matrix opC(igrid.size(), nN);
      wnl.compute<Internal::MHDFloat>(opC, nN, l, igrid, iweights.array()*igrid.array().pow(-1), ev::Set());

      tOp = (opC.transpose()*opB*opA.transpose()).transpose();
#endif
      op = tOp.cast<MHDFloat>().leftCols(nPoly);

      assert(op.rows() == igrid.size());
      assert(op.cols() == nPoly);
   }

   void DivR1D1R1::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}

#undef QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
