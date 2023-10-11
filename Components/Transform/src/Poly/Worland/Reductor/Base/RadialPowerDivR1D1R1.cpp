/**
 * @file RadialPowerDivR1D1R1.cpp
 * @brief Source of the implementation of the Worland 1/R D R radial power spectrum operator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/RadialPowerDivR1D1R1.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   RadialPowerDivR1D1R1<base_t>::RadialPowerDivR1D1R1()
      : IWorlandRadialPower()
   {
      this->setProfileTag();
   }

   void RadialPowerDivR1D1R1<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      if(l == 0)
      {
         op.setZero();
      }
      else
      {
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::r_1drWnl wnl;
         wnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());
      }
   }

   void RadialPowerDivR1D1R1<base_t>::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_REDUIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_REDUIMPL_OTF
         int l = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fastSize(i);
         if(l == 0)
         {
            rOut.setZero();
         }
         else
         {
            namespace ev = Polynomial::Worland::Evaluator;
            Polynomial::Worland::r_1drWnl wnl;
            MatrixZ tmp(this->mGrid.size(), in.cols());
            wnl.compute<MHDComplex>(tmp, nPoly, l, this->mGrid, Internal::Array(), ev::OuterProduct(in));
            rOut = tmp.array().abs2();
         }
      #endif //defined QUICC_WORLAND_REDUIMPL_MATRIX
   }

} // Reductor
} // Worland
} // Poly
} // Transform
} // QuICC
