/**
 * @file PowerSLAPLR2.cpp
 * @brief Source of the implementation of the Worland Spherical Laplacian R^2 power spectrum operator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Reductor/PowerSLAPLR2.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/slaplWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   PowerSLAPLR2::PowerSLAPLR2()
      : IWorlandPower(1)
   {
      this->setProfileTag();
   }

   PowerSLAPLR2::~PowerSLAPLR2()
   {
   }

   void PowerSLAPLR2::makeOperator(Matrix& op, Matrix& eop, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fastSize(i);

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::slaplWnl bwnl;
      bwnl.compute<MHDFloat>(op, nPoly, l, igrid, internal::Array(), ev::Set());

      Polynomial::Worland::Wnl fwnl(Polynomial::Worland::Wnl::ALPHA_SPHENERGY,Polynomial::Worland::Wnl::DBETA_SPHENERGY);

      eop.resize(igrid.size(), nPoly);
      fwnl.compute<MHDFloat>(eop, nPoly, l, igrid, iweights, ev::Set());
   }

   void PowerSLAPLR2::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_REDUIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_REDUIMPL_OTF
         int l = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fastSize(i);

         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::slaplWnl bwnl;
         MatrixZ tmp(this->mGrid.size(), in.cols());
         bwnl.compute<MHDComplex>(tmp, nPoly, l, this->mGrid, internal::Array(), ev::OuterProduct<MHDComplex>(in));

         Polynomial::Worland::Wnl fwnl(Polynomial::Worland::Wnl::ALPHA_SPHENERGY,Polynomial::Worland::Wnl::DBETA_SPHENERGY);
         MatrixZ tmpB(nPoly, in.cols());
         fwnl.compute<MHDComplex>(tmpB, nPoly, l, this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(tmp));
         rOut = tmpB.array().abs2();
      #endif //defined QUICC_WORLAND_REDUIMPL_MATRIX
   }

}
}
}
}
}
