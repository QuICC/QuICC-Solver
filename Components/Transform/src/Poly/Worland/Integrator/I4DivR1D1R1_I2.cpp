/**
 * @file I4DivR1D1R1_I2.cpp
 * @brief Source of the implementation of the Worland 1/R1 D R1 integrator but 0 mode is I2 P integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Integrator/I4DivR1D1R1_I2.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/dWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   I4DivR1D1R1_I2::I4DivR1D1R1_I2()
   {
   }

   I4DivR1D1R1_I2::~I4DivR1D1R1_I2()
   {
   }

   void I4DivR1D1R1_I2::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      if(l == 0)
      {
         l = 1;
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;

         op.resize(igrid.size(), nPoly);
         wnl.compute<MHDFloat>(op, nPoly, l, igrid, iweights, ev::Set());

         MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
         MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
         ::QuICC::SparseSM::Worland::I2 spasm(nPoly, nPoly, a, b, l);
         op = (spasm.mat()*op.transpose()).transpose();
      } else
      {
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;
#if 0
         // **************************************************
         // Formulation without explicit grid:
         // Operates on polynomials with l = l-1
         int l_in = std::abs(l-1);
         int n_in = nPoly + 1;

         Matrix opA(igrid.size(), n_in);
         wnl.compute<MHDFloat>(opA, n_in, l_in, igrid, iweights, ev::Set());

         Matrix opB(igrid.size(), n_in);
         Polynomial::Worland::r_1drWnl r_1drWnl;
         r_1drWnl.compute<MHDFloat>(opB, n_in, l_in, igrid, internal::Array(), ev::Set());

         Matrix opC(igrid.size(), nPoly);
         Polynomial::Worland::Wnl wnlB;
         wnlB.compute<MHDFloat>(opC, nPoly, l, igrid, iweights, ev::Set());

         op = (opC.transpose()*opB*opA.transpose()).transpose();
#else

         // **************************************************
         // Alternative formulation of operators:
         // This version uses explicit radial factor to work on l polynomials

         Matrix opA(igrid.size(), nPoly);
         wnl.compute<MHDFloat>(opA, nPoly, l, igrid, iweights, ev::Set());

         Matrix opB(igrid.size(), nPoly);
         Polynomial::Worland::dWnl dWnl;
         dWnl.compute<MHDFloat>(opB, nPoly, l, igrid, internal::Array(), ev::Set());

         op = (opA.transpose()*igrid.cast<MHDFloat>().array().pow(-1).matrix().asDiagonal()*opB*opA.transpose()*igrid.cast<MHDFloat>().asDiagonal()).transpose();
#endif

         // Multiply by Quasi-inverse
         MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
         MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
         ::QuICC::SparseSM::Worland::I4 spasm(nPoly, nPoly, a, b, l);
         op = (spasm.mat()*op.transpose()).transpose();

         assert(op.rows() == igrid.size());
         assert(op.cols() == nPoly);
      }
   }

   void I4DivR1D1R1_I2::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_INTGIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_INTGIMPL_OTF
         int l = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fastSize(i);

         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;
         if(l == 0)
         {
            l = 1;
            wnl.compute<MHDComplex>(rOut, nPoly, l, this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(in));

            MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
            MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
            ::QuICC::SparseSM::Worland::I2 spasm(nPoly, nPoly, a, b, l);
            rOut = spasm.mat()*rOut;
         } else
         {
            int l_in = std::abs(l-1);

            wnl.compute<MHDComplex>(rOut, nPoly, l_in, this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(in));

            MatrixZ tmp(in.rows(), in.cols());
            Polynomial::Worland::r_1drWnl r_1drWnl;
            r_1drWnl.compute<MHDComplex>(tmp, nPoly, l_in, this->mGrid, internal::Array(), ev::OuterProduct<MHDComplex>(rOut));

            wnl.compute<MHDComplex>(rOut, nPoly, l, this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(tmp));

            MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
            MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
            ::QuICC::SparseSM::Worland::I4 spasm(nPoly, nPoly, a, b, l);
            rOut = spasm.mat()*rOut;
         }
      #endif //defined QUICC_WORLAND_INTGIMPL_MATRIX
   }

}
}
}
}
}
