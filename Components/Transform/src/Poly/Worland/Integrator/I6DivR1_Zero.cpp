/**
 * @file I6DivR1_Zero.cpp
 * @brief Source of the implementation of the Worland 1/R1 integrator but 0 mode is zeroed
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Integrator/I6DivR1_Zero.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"
#include "QuICC/SparseSM/Worland/I6.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   I6DivR1_Zero::I6DivR1_Zero()
   {
   }

   I6DivR1_Zero::~I6DivR1_Zero()
   {
   }

   void I6DivR1_Zero::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      if(l == 0)
      {
         op.setZero();
      } else
      {
         Polynomial::Worland::Wnl wnl;

         wnl.compute<MHDFloat>(op, nPoly, std::abs(l-1), igrid, iweights, ev::Set());

         Matrix opA(igrid.size(), nPoly);
         Polynomial::Worland::r_1Wnl r_1Wnl;
         r_1Wnl.compute<MHDFloat>(opA, nPoly, std::abs(l-1), igrid, internal::Array(), ev::Set());

         Matrix opB(igrid.size(), nPoly);
         Polynomial::Worland::Wnl wnlB;
         wnlB.compute<MHDFloat>(opB, nPoly, l, igrid, iweights, ev::Set());

         MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
         MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
         ::QuICC::SparseSM::Worland::I6 spasm(nPoly, nPoly, a, b, l);
         op = (spasm.mat()*opB.transpose()*opA*op.transpose()).transpose();
      }
   }

   void I6DivR1_Zero::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
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
            namespace ev = Polynomial::Worland::Evaluator;
            Polynomial::Worland::Wnl wnl;
            wnl.compute<MHDComplex>(rOut, nPoly, std::abs(l-1), this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(in));
            MatrixZ tmp(in.rows(), in.cols());
            Polynomial::Worland::r_1Wnl r_1Wnl;
            r_1Wnl.compute<MHDComplex>(tmp, nPoly, std::abs(l-1), this->mGrid, internal::Array(), ev::OuterProduct<MHDComplex>(rOut));
            wnl.compute<MHDComplex>(rOut, nPoly, l, this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(tmp));

            MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
            MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
            ::QuICC::SparseSM::Worland::I6 spasm(nPoly, nPoly, a, b, l);
            rOut = spasm.mat()*rOut;
         }
      #endif //defined QUICC_WORLAND_INTGIMPL_MATRIX
   }

}
}
}
}
}
