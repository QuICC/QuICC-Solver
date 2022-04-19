/**
 * @file I6CylLaplh_I4D1R1.cpp
 * @brief Source of the implementation of the Worland 1/R1 D R1 integrator but 0 mode is I4 D R1 integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Integrator/I6CylLaplh_I4D1R1.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/SparseSM/Worland/I4D1R1.hpp"
#include "QuICC/SparseSM/Worland/I6CylLaplh.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   I6CylLaplh_I4D1R1::I6CylLaplh_I4D1R1()
   {
   }

   I6CylLaplh_I4D1R1::~I6CylLaplh_I4D1R1()
   {
   }

   void I6CylLaplh_I4D1R1::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;
      wnl.compute<MHDFloat>(op, nPoly, l, igrid, iweights, ev::Set());
      if(l == 0)
      {
         MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
         MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
         ::QuICC::SparseSM::Worland::I4D1R1 spasm(nPoly, nPoly, a, b, 1);
         op = (spasm.mat()*op.transpose()).transpose();
      } else
      {
         MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
         MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
         ::QuICC::SparseSM::Worland::I6CylLaplh spasm(nPoly, nPoly, a, b, l);
         op = (spasm.mat()*op.transpose()).transpose();
      }
   }

   void I6CylLaplh_I4D1R1::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_INTGIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_INTGIMPL_OTF
         int l = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fastSize(i);
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;
         wnl.compute<MHDComplex>(rOut, nPoly,l, this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(in));
         if(l == 0)
         {
            MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
            MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
            ::QuICC::SparseSM::Worland::I4D1R1 spasm(nPoly, nPoly, a, b, 1);
            rOut = spasm.mat()*rOut;
         } else
         {
            MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
            MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
            ::QuICC::SparseSM::Worland::I6CylLaplh spasm(nPoly, nPoly, a, b, l);
            rOut = spasm.mat()*rOut;
         }
      #endif //defined QUICC_WORLAND_INTGIMPL_MATRIX
   }

}
}
}
}
}
