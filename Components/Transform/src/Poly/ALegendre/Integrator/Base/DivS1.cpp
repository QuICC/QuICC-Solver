/**
 * @file DivS1.cpp
 * @brief Source of the implementation of the associated Legendre P integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/DivS1.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/InnerProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   void DivS1<base_t>::makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const
   {
      int m = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::ALegendre::Evaluator;
      Polynomial::ALegendre::sin_1Plm dplm;
      dplm.compute<MHDFloat>(op, nPoly, m, igrid, iweights, ev::Set());
   }

   void DivS1<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {

      #if defined QUICC_ALEGENDRE_INTGIMPL_OTF
         int m = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;
         namespace ev = Polynomial::ALegendre::Evaluator;
         Polynomial::ALegendre::sin_1Plm dplm;
         dplm.compute<MHDComplex>(rOut, nPoly, m, this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(in));
      #elif defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
         rOut = this->mOps.at(i).transpose()*in;
      #endif //defined QUICC_ALEGENDRE_INTGIMPL_OTF

   }

}
}
}
}
}
