/**
 * @file DivS1.cpp
 * @brief Source of the implementation of the associated Legendre 1/Sin projector
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Projector/DivS1.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   DivS1::DivS1()
   {
   }

   DivS1::~DivS1()
   {
   }

   void DivS1::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int m = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::ALegendre::Evaluator;
      Polynomial::ALegendre::sin_1Plm plm;
      plm.compute<MHDFloat>(op, nPoly, m, igrid, internal::Array(), ev::Set());
   }

   void DivS1::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_ALEGENDRE_PROJIMPL_OTF
         int m = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;
         namespace ev = Polynomial::ALegendre::Evaluator;
         Polynomial::ALegendre::sin_1Plm plm;
         plm.compute<MHDComplex>(rOut, nPoly, m, this->mGrid, internal::Array(), ev::OuterProduct<MHDComplex>(in));
      #elif defined QUICC_ALEGENDRE_PROJIMPL_MATRIX
         rOut = this->mOps.at(i).transpose()*in;
      #endif //defined QUICC_ALEGENDRE_PROJIMPL_OTF
   }

}
}
}
}
}
