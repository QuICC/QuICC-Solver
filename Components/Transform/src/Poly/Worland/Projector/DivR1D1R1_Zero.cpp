/**
 * @file DivR1D1R1_Zero.cpp
 * @brief Source of the implementation of the Worland 1/R D R projector and zero for l = 0
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/BreakPoint.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/DivR1D1R1_Zero.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   DivR1D1R1_Zero::DivR1D1R1_Zero()
   {
      this->mProfileId = Debug::Profiler::WORLANDPROJ_DIVR1D1R1_ZERO;
   }

   DivR1D1R1_Zero::~DivR1D1R1_Zero()
   {
   }

   void DivR1D1R1_Zero::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
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
         Polynomial::Worland::r_1drWnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, l, igrid, internal::Array(), ev::Set());
      }
   }

   void DivR1D1R1_Zero::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
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
            Polynomial::Worland::r_1drWnl wnl;
            wnl.compute<MHDComplex>(rOut, this->mspSetup->fastSize(i), this->mspSetup->slow(i), this->mGrid, internal::Array(), ev::OuterProduct(in));
         }
      #endif //defined QUICC_WORLAND_PROJIMPL_MATRIX
   }

}
}
}
}
}
