/**
 * @file RadialPower.cpp
 * @brief Source of the implementation of the Worland power spectrum operator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/RadialPower.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   RadialPower<base_t>::RadialPower()
      : IWorlandRadialPower()
   {
      this->setProfileTag();
   }

   void RadialPower<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;
      wnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());
   }

   void RadialPower<base_t>::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

} // Reductor
} // Worland
} // Poly
} // Transform
} // QuICC
