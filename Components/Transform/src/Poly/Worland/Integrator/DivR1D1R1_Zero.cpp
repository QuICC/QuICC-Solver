/**
 * @file DivR1D1R1_Zero.cpp
 * @brief Source of the implementation of the Worland 1/R D R integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Integrator/DivR1D1R1_Zero.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/dWnl.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   DivR1D1R1_Zero::DivR1D1R1_Zero()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_DIVR1D1R1;
   }

   DivR1D1R1_Zero::~DivR1D1R1_Zero()
   {
   }

   void DivR1D1R1_Zero::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
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

         internal::Matrix tOp(igrid.size(), nN);
#if 1
         // **************************************************
         // Formulation without explicit grid:
         // Operates on polynomials with l = l-1
         int l_in = std::abs(l-1);
         int l_out = std::abs(l-0);
         int n_in = nN + 1;
         this->checkGridSize(n_in, l_in, igrid.size());

         internal::Matrix opA(igrid.size(), n_in);
         wnl.compute<internal::MHDFloat>(opA, n_in, l_in, igrid, iweights, ev::Set());

         internal::Matrix opB(igrid.size(), n_in);
         Polynomial::Worland::r_1drWnl r_1drWnl;
         r_1drWnl.compute<internal::MHDFloat>(opB, n_in, l_in, igrid, internal::Array(), ev::Set());

         internal::Matrix opC(igrid.size(), nN);
         Polynomial::Worland::Wnl wnlB;
         wnlB.compute<internal::MHDFloat>(opC, nN, l_out, igrid, iweights, ev::Set());

         tOp = (opC.transpose()*opB*opA.transpose()).transpose();
#else
         // **************************************************
         // Alternative formulation of operators:
         // This version uses explicit radial factor to work on l polynomials

         internal::Matrix opA(igrid.size(), nN);
         wnl.compute<internal::MHDFloat>(opA, nN, l, igrid, iweights.array()*igrid.array(), ev::Set());

         internal::Matrix opB(igrid.size(), nN);
         Polynomial::Worland::dWnl dWnl;
         dWnl.compute<internal::MHDFloat>(opB, nN, l, igrid, internal::Array(), ev::Set());

         internal::Matrix opC(igrid.size(), nN);
         wnl.compute<internal::MHDFloat>(opC, nN, l, igrid, iweights.array()*igrid.array().pow(-1), ev::Set());

         tOp = (opC.transpose()*opB*opA.transpose()).transpose();
#endif
         op = tOp.cast<MHDFloat>().leftCols(nPoly);

         assert(op.rows() == igrid.size());
         assert(op.cols() == nPoly);
      }
   }

   void DivR1D1R1_Zero::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_INTGIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_INTGIMPL_OTF
         int l = this->mspSetup->slow(i);
         if(l == 0)
         {
            rOut.setZero();
         } else
         {
            int nPoly = this->mspSetup->fastSize(i);
            int l_in = std::abs(l-1);

            namespace ev = Polynomial::Worland::Evaluator;
            Polynomial::Worland::Wnl wnl;

            wnl.compute<MHDComplex>(rOut, nPoly, l_in, this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(in));

            MatrixZ tmp(in.rows(), in.cols());
            Polynomial::Worland::r_1drWnl r_1drWnl;
            r_1drWnl.compute<MHDComplex>(tmp, nPoly, l_in, this->mGrid, internal::Array(), ev::OuterProduct<MHDComplex>(rOut));

            wnl.compute<MHDComplex>(rOut, nPoly, l, this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(tmp));
         }
      #endif //defined QUICC_WORLAND_INTGIMPL_MATRIX
   }

}
}
}
}
}
