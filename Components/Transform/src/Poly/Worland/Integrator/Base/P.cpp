/**
 * @file P.cpp
 * @brief Source of the implementation of the Worland P integrator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/P.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   P<base_t>::P()
      : IWorlandIntegrator()
   {
      this->setProfileTag();
   }

   void P<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;

      // Internal computation uses dealiased modes
      int nN = nPoly;
      this->checkGridSize(nN, l, igrid.size());

      Internal::Matrix tOp(igrid.size(), nN);

      wnl.compute<Internal::MHDFloat>(tOp, nN, l, igrid, iweights, ev::Set());

      op = tOp.cast<MHDFloat>().leftCols(nPoly);
   }

   void P<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_INTGIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_INTGIMPL_OTF
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;
         wnl.compute<MHDComplex>(rOut, this->mspSetup->fastSize(i), this->mspSetup->slow(i), this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(in));
      #endif //defined QUICC_WORLAND_INTGIMPL_MATRIX
   }

}
}
}
}
}
