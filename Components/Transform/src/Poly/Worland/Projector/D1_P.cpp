/**
 * @file D1_P.cpp
 * @brief Source of the implementation of the Worland D projector but 0 mode is P projector with l = 1
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Projector/D1_P.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/dWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   D1_P::D1_P()
   {
   }

   D1_P::~D1_P()
   {
   }

   void D1_P::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      if(l == 0)
      {
         Polynomial::Worland::Wnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, 1, igrid, internal::Array(), ev::Set());
      } else
      {
         Polynomial::Worland::dWnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, l, igrid, internal::Array(), ev::Set());
      }
   }

   void D1_P::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_PROJIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_PROJIMPL_OTF
         int l = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fastSize(i);
         namespace ev = Polynomial::Worland::Evaluator;
         if(l == 0)
         {
            Polynomial::Worland::Wnl wnl;
            wnl.compute<MHDComplex>(rOut, nPoly, 1, this->mGrid, internal::Array(), ev::OuterProduct(in));
         } else
         {
            Polynomial::Worland::dWnl wnl;
            wnl.compute<MHDComplex>(rOut, nPoly, l, this->mGrid, internal::Array(), ev::OuterProduct(in));
         }
      #endif //defined QUICC_WORLAND_PROJIMPL_MATRIX
   }

}
}
}
}
}
