/**
 * @file P.cpp
 * @brief Source of the implementation of the associated Legendre P integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/P.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/InnerProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   P::P()
   {
   }

   P::~P()
   {
   }

   void P::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int m = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::ALegendre::Evaluator;
      Polynomial::ALegendre::Plm plm;
      plm.compute<MHDFloat>(op, nPoly, m, igrid, iweights, ev::Set());
   }

   void P::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_ALEGENDRE_INTGIMPL_OTF
         int m = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1;
         namespace ev = Polynomial::ALegendre::Evaluator;
         Polynomial::ALegendre::Plm plm;
         plm.compute<MHDComplex>(rOut, nPoly, m, this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(in));
      #elif defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
         rOut = this->mOps.at(i).transpose()*in;
      #endif //defined QUICC_ALEGENDRE_INTGIMPL_OTF
   }

}
}
}
}
}
