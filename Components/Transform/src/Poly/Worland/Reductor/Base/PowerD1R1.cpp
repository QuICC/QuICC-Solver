/**
 * @file PowerD1R1.cpp
 * @brief Source of the implementation of the Worland D R power spectrum operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/PowerD1R1.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   PowerD1R1<base_t>::PowerD1R1()
      : IWorlandPower(0)
   {
      this->setProfileTag();
   }

   void PowerD1R1<base_t>::makeOperator(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fastSize(i);

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::r_1drWnl bwnl;
      bwnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());

      eop.resize(igrid.size(), nPoly);
      Polynomial::Worland::Wnl fwnl(Polynomial::Worland::Wnl::ALPHA_SPHENERGY,Polynomial::Worland::Wnl::DBETA_SPHENERGY);
      if(l == 0)
      {
         eop.setZero();
      }
      else
      {
         fwnl.compute<MHDFloat>(eop, nPoly, std::abs(l-1), igrid, iweights, ev::Set());
      }
   }

   void PowerD1R1<base_t>::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_REDUIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_REDUIMPL_OTF
         int l = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fastSize(i);

         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::r_1drWnl bwnl;
         MatrixZ tmp(this->mGrid.size(), in.cols());
         bwnl.compute<MHDComplex>(tmp, nPoly, l, this->mGrid, Internal::Array(), ev::OuterProduct<MHDComplex>(in));

         Polynomial::Worland::Wnl fwnl(Polynomial::Worland::Wnl::ALPHA_SPHENERGY,Polynomial::Worland::Wnl::DBETA_SPHENERGY);
         MatrixZ tmpB(nPoly, in.cols());
         fwnl.compute<MHDComplex>(tmpB, nPoly, l-1, this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(tmp));
         rOut = tmpB.array().abs2();
      #endif //defined QUICC_WORLAND_REDUIMPL_MATRIX
   }

}
}
}
}
}
