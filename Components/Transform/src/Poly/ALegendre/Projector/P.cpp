/**
 * @file P.cpp
 * @brief Source of the implementation of the associated Legendre P projector
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Projector/P.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template<typename OpTypes>
   P<OpTypes>::P()
   {
   }

   template<typename OpTypes>
   P<OpTypes>::~P()
   {
   }

   template<typename OpTypes>
   void P<OpTypes>::makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const
   {
      int m = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::ALegendre::Evaluator;
      Polynomial::ALegendre::Plm plm;
      plm.compute<MHDFloat>(op, nPoly, m, igrid, OpArray(), ev::Set());
   }

   template<typename OpTypes>
   void P<OpTypes>::applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const
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

   template class P<IALegendreOperatorTypes>;
   template class P<PIALegendreOperatorTypes>;
}
}
}
}
}
