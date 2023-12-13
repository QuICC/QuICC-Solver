/**
 * @file DivR1D1R1_Zero.cpp
 * @brief Source of the implementation of the Worland 1/R D R integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/DivR1D1R1_Zero.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

// It is not clear yet which implementation is more accurate
#define QUICC_AVOID_EXPLICIT_RADIAL_FACTOR

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   DivR1D1R1_Zero<base_t>::DivR1D1R1_Zero()
   {
      this->setProfileTag();
   }

   void DivR1D1R1_Zero<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
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
         // Internal computation uses dealiased modes
         this->checkGridSize(nPoly, l, igrid.size());
#ifdef QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
         using poly_t = QuICC::Polynomial::Worland::implicit_t;
#else
         using poly_t = QuICC::Polynomial::Worland::explicit_t;
#endif
         Polynomial::Worland::r_1drWnl<poly_t> wnl;
         wnl.compute<MHDFloat>(op, nPoly, l, igrid, iweights, ev::Set());
      }
   }

   void DivR1D1R1_Zero<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
#undef QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
