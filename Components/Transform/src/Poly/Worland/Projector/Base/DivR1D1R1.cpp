/**
 * @file DivR1D1R1.cpp
 * @brief Source of the implementation of the Worland 1/R D R projector
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/DivR1D1R1.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   DivR1D1R1<base_t>::DivR1D1R1()
   {
      this->setProfileTag();
   }

   void DivR1D1R1<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::r_1drWnl wnl;
      wnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());
   }

   void DivR1D1R1<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_PROJIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_PROJIMPL_OTF
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::r_1drWnl wnl;
         wnl.compute<MHDComplex>(rOut, this->mspSetup->fastSize(i), this->mspSetup->slow(i), this->mGrid, Internal::Array(), ev::OuterProduct(in));
      #endif //defined QUICC_WORLAND_PROJIMPL_MATRIX
   }

}
}
}
}
}
