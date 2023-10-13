/**
 * @file CylLaplh_DivR1D1R1.cpp
 * @brief Source of the implementation of the Worland cylindrical horizontal laplacian projector but 0 mode is 1/R D R projector
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/CylLaplh_DivR1D1R1.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/claplhWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   void CylLaplh_DivR1D1R1<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      if(l == 0)
      {
         Polynomial::Worland::r_1drWnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, 1, igrid, Internal::Array(), ev::Set());
      } else
      {
         Polynomial::Worland::claplhWnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());
      }
   }

   void CylLaplh_DivR1D1R1<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
