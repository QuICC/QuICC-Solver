/**
 * @file D1CylLaplh_D1DivR1D1R1.cpp
 * @brief Source of the implementation of the Worland D of cylindrical horizontal laplacian projector but 0 mode is D 1/R D R projector with l = 1
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Projector/D1CylLaplh_D1DivR1D1R1.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/dr_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/dclaplhWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   D1CylLaplh_D1DivR1D1R1::D1CylLaplh_D1DivR1D1R1()
   {
   }

   D1CylLaplh_D1DivR1D1R1::~D1CylLaplh_D1DivR1D1R1()
   {
   }

   void D1CylLaplh_D1DivR1D1R1::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      if(l == 0)
      {
         l = 1;
         Polynomial::Worland::dr_1drWnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, l, igrid, internal::Array(), ev::Set());
      } else
      {
         Polynomial::Worland::dclaplhWnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, l, igrid, internal::Array(), ev::Set());
      }
   }

   void D1CylLaplh_D1DivR1D1R1::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_PROJIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_PROJIMPL_OTF
         int l = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fastSize(i);
         namespace ev = Polynomial::Worland::Evaluator;
         if(l == 0)
         {
            l = 1;
            Polynomial::Worland::dr_1drWnl wnl;
            wnl.compute<MHDComplex>(rOut, nPoly, l, this->mGrid, internal::Array(), ev::OuterProduct(in));
         } else
         {
            Polynomial::Worland::dclaplhWnl wnl;
            wnl.compute<MHDComplex>(rOut, nPoly, l, this->mGrid, internal::Array(), ev::OuterProduct(in));
         }
      #endif //defined QUICC_WORLAND_PROJIMPL_MATRIX
   }

}
}
}
}
}
