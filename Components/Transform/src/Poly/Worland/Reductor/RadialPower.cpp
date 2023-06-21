/**
 * @file RadialPower.cpp
 * @brief Source of the implementation of the Worland power spectrum operator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/RadialPower.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   RadialPower::RadialPower()
      : IWorlandRadialPower()
   {
      this->setProfileTag();
   }

   void RadialPower::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;
      wnl.compute<MHDFloat>(op, nPoly, l, igrid, internal::Array(), ev::Set());
   }

   void RadialPower::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_REDUIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_REDUIMPL_OTF
         int l = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fastSize(i);

         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;
         MatrixZ tmp(this->mGrid.size(), in.cols());
         wnl.compute<MHDComplex>(tmp, nPoly, l, this->mGrid, internal::Array(), ev::OuterProduct(in));
         rOut = tmp.array().abs2();
      #endif //defined QUICC_WORLAND_REDUIMPL_MATRIX
   }

} // Reductor
} // Worland
} // Poly
} // Transform
} // QuICC
