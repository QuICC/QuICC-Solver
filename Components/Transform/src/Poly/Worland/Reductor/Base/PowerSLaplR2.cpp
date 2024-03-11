/**
 * @file PowerSLaplR2.cpp
 * @brief Source of the implementation of the Worland Spherical Laplacian R^2 power spectrum operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/PowerSLaplR2.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
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

   PowerSLaplR2<base_t>::PowerSLaplR2()
      : IWorlandPower(1)
   {
      this->setProfileTag();
   }

   void PowerSLaplR2<base_t>::makeOperator(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fastSize(i);

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::slaplWnl bwnl;
      bwnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());

      Polynomial::Worland::Wnl fwnl(Polynomial::Worland::worland_sphenergy_t::ALPHA,Polynomial::Worland::worland_sphenergy_t::DBETA);

      eop.resize(igrid.size(), nPoly);
      fwnl.compute<MHDFloat>(eop, nPoly, l, igrid, iweights, ev::Set());
   }

   void PowerSLaplR2<base_t>::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
