/**
 * @file I2_Zero.cpp
 * @brief Source of the implementation of the Worland P integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Integrator/I2_Zero.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   I2_Zero::I2_Zero()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_I2_ZERO;
   }

   I2_Zero::~I2_Zero()
   {
   }

   void I2_Zero::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      if(l == 0)
      {
         op.setZero();
      } else
      {
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, l, igrid, iweights, ev::Set());

         MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
         MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
         ::QuICC::SparseSM::Worland::I2 spasm(nPoly, nPoly, a, b, l);
         op = (spasm.mat()*op.transpose()).transpose();
      }
   }

   void I2_Zero::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
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
            wnl.compute<MHDComplex>(rOut, nPoly, this->mspSetup->slow(i), this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(in));

            MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
            MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
            ::QuICC::SparseSM::Worland::I2 spasm(nPoly, nPoly, a, b, l);
            rOut = spasm.mat()*rOut;
         }
      #endif //defined QUICC_WORLAND_INTGIMPL_MATRIX
   }

}
}
}
}
}
