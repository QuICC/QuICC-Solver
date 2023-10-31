/**
 * @file Power.cpp
 * @brief Source of the implementation of the Worland power spectrum operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/Power.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   Power<base_t>::Power()
      : IWorlandPower(1)
   {
      this->setProfileTag();
   }

   void Power<base_t>::makeOperator(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fastSize(i);

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::r_1Wnl<QuICC::Polynomial::Worland::recurrence_t> bwnl;
      bwnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());

      eop.resize(igrid.size(), nPoly);
      Polynomial::Worland::Wnl fwnl(Polynomial::Worland::Wnl::ALPHA_SPHENERGY,Polynomial::Worland::Wnl::DBETA_SPHENERGY);
      if(l == 0)
      {
         eop.setZero();
      }
      else
      {
         fwnl.compute<MHDFloat>(eop, nPoly, l-1, igrid, iweights, ev::Set());
      }
   }

   void Power<base_t>::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
