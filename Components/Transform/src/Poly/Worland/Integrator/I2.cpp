/**
 * @file I2.cpp
 * @brief Source of the implementation of the Worland I2 integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Integrator/I2.hpp"

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

   I2::I2()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_I2;
   }

   I2::~I2()
   {
   }

   void I2::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);

      // Internal computation uses dealiased modes
      const int extraN = 3; // I2 has 3 superdiagonals
      int nN = nPoly + extraN;
      this->checkGridSize(nN, l, igrid.size());

      internal::Matrix tOp(igrid.size(), nN);

      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;
      wnl.compute<internal::MHDFloat>(tOp, nN, l, igrid, iweights, ev::Set());

      auto a = wnl.alpha(l);
      auto b = wnl.dBeta();
      ::QuICC::SparseSM::Worland::I2 spasm(nN, nN, a, b, l);
      tOp = (spasm.mat()*tOp.transpose()).transpose();
      op = tOp.cast<MHDFloat>().leftCols(nPoly);
   }

   void I2::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_INTGIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_INTGIMPL_OTF
         int l = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fastSize(i);
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;
         wnl.compute<MHDComplex>(rOut, nPoly, this->mspSetup->slow(i), this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(in));

         MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
         MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
         ::QuICC::SparseSM::Worland::I2 spasm(nPoly, nPoly, a, b, l);
         rOut = spasm.mat()*rOut;
      #endif //defined QUICC_WORLAND_INTGIMPL_MATRIX
   }

}
}
}
}
}
