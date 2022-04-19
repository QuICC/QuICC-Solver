/**
 * @file DivR1CylLaplh_Zero.cpp
 * @brief Source of the implementation of the Worland D of horizontal cylindrical laplacian projector but 0 mode is zeroed
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Projector/DivR1CylLaplh_Zero.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/r_1claplhWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   DivR1CylLaplh_Zero::DivR1CylLaplh_Zero()
   {
   }

   DivR1CylLaplh_Zero::~DivR1CylLaplh_Zero()
   {
   }

   void DivR1CylLaplh_Zero::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      if(l == 0)
      {
         op.setZero();
      } else
      {
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::r_1claplhWnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, l, igrid, internal::Array(), ev::Set());
      }
   }

   void DivR1CylLaplh_Zero::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_PROJIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_PROJIMPL_OTF
         int l = this->mspSetup->slow(i);
         if(l == 0)
         {
            rOut.setZero();
         } else
         {
            namespace ev = Polynomial::Worland::Evaluator;
            Polynomial::Worland::r_1claplhWnl wnl;
            wnl.compute<MHDComplex>(rOut, this->mspSetup->fastSize(i), l, this->mGrid, internal::Array(), ev::OuterProduct(in));
         }
      #endif //defined QUICC_WORLAND_PROJIMPL_MATRIX
   }

}
}
}
}
}
