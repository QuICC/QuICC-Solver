/**
 * @file CylLaplh_DivR1D1R1.cpp
 * @brief Source of the implementation of the Worland cylindrical horizontal laplacian projector but 0 mode is 1/R D R projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Projector/CylLaplh_DivR1D1R1.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/claplhWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   CylLaplh_DivR1D1R1::CylLaplh_DivR1D1R1()
   {
   }

   CylLaplh_DivR1D1R1::~CylLaplh_DivR1D1R1()
   {
   }

   void CylLaplh_DivR1D1R1::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      if(l == 0)
      {
         Polynomial::Worland::r_1drWnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, 1, igrid, internal::Array(), ev::Set());
      } else
      {
         Polynomial::Worland::claplhWnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, l, igrid, internal::Array(), ev::Set());
      }
   }

   void CylLaplh_DivR1D1R1::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_PROJIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_PROJIMPL_OTF
         int l = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fastSize(i);
         namespace ev = Polynomial::Worland::Evaluator;
         if(l == 0)
         {
            Polynomial::Worland::r_1drWnl wnl;
            wnl.compute<MHDComplex>(rOut, nPoly, 1, this->mGrid, internal::Array(), ev::OuterProduct(in));
         } else
         {
            Polynomial::Worland::claplhWnl wnl;
            wnl.compute<MHDComplex>(rOut, nPoly, l, this->mGrid, internal::Array(), ev::OuterProduct(in));
         }
      #endif //defined QUICC_WORLANd_PROJIMPL_MATRIX
   }

}
}
}
}
}
