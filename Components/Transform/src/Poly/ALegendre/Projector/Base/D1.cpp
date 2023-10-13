/**
 * @file D1.cpp
 * @brief Source of the implementation of the associated Legendre D projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/D1.hpp"
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   void D1<base_t>::makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const
   {
      int m = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::ALegendre::Evaluator;
      Polynomial::ALegendre::dPlm dplm;
      dplm.compute<MHDFloat>(op, nPoly, m, igrid, OpArray(), ev::Set());
   }

   void D1<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      rOut = this->mOps.at(i).transpose()*in;
   }

}
}
}
}
}
