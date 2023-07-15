/**
 * @file P.cpp
 * @brief Source of the implementation of the associated Legendre P projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/P.hpp"
#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   void P<base_t>::makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const
   {
      int m = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::ALegendre::Evaluator;
      Polynomial::ALegendre::Plm plm;
      plm.compute<MHDFloat>(op, nPoly, m, igrid, OpArray(), ev::Set());
   }

   void P<base_t>::applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const
   {
      #if defined QUICC_ALEGENDRE_PROJIMPL_OTF
         int m = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;
         namespace ev = Polynomial::ALegendre::Evaluator;
         Polynomial::ALegendre::Plm plm;
         plm.compute<MHDComplex>(rOut, nPoly, m, this->mGrid, OpArray(), ev::OuterProduct<MHDComplex>(in));
      #elif defined QUICC_ALEGENDRE_PROJIMPL_MATRIX
         rOut = this->mOps.at(i).transpose()*in;
      #endif //defined QUICC_ALEGENDRE_PROJIMPL_OTF
   }

}
}
}
}
}
