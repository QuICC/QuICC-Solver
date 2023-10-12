/**
 * @file D1_P.cpp
 * @brief Source of the implementation of the Worland D projector but 0 mode is P projector with l = 1
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/D1_P.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/dWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   void D1_P<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      if(l == 0)
      {
         Polynomial::Worland::Wnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, 1, igrid, Internal::Array(), ev::Set());
      } else
      {
         Polynomial::Worland::dWnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());
      }
   }

   void D1_P<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
