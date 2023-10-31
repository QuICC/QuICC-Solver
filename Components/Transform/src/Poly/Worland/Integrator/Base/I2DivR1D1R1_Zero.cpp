/**
 * @file I2DivR1D1R1_Zero.cpp
 * @brief Source of the implementation of the Worland I2 1/R D R integrator but 0 mode is zeroed
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/I2DivR1D1R1_Zero.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/dWnl.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"

// It is not clear yet which implementation is more accurate
#define QUICC_AVOID_EXPLICIT_RADIAL_FACTOR

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   I2DivR1D1R1_Zero<base_t>::I2DivR1D1R1_Zero()
   {
      this->setProfileTag();
   }

   void I2DivR1D1R1_Zero<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      if(l == 0)
      {
         op.resize(igrid.size(), nPoly);
         op.setZero();
      } else
      {
         // Internal computation uses dealiased modes
         const int extraN = 3*(!this->mcTruncQI); // I2 has 3 superdiagonals
         int nN = nPoly + extraN;
         this->checkGridSize(nN, l, igrid.size());

         Internal::Matrix tOp(igrid.size(), nN);

         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;

#ifdef QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
         using poly_t = QuICC::Polynomial::Worland::implicit_t;
#else
         using poly_t = QuICC::Polynomial::Worland::explicit_t;
#endif
         Polynomial::Worland::r_1drWnl<poly_t> r_1drWnl;
         r_1drWnl.compute<Internal::MHDFloat>(tOp, nN, l, igrid, iweights, ev::Set());

         // Multiply by Quasi-inverse
         auto a = wnl.alpha(l);
         auto b = wnl.dBeta();
         ::QuICC::SparseSM::Worland::I2 spasm(nN, nN, a, b, l, 1*this->mcTruncQI);
         tOp = (spasm.mat()*tOp.transpose()).transpose();
         op = tOp.cast<MHDFloat>().leftCols(nPoly);

         assert(op.rows() == igrid.size());
         assert(op.cols() == nPoly);
      }
   }

   void I2DivR1D1R1_Zero<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}

#undef QUICC_AVOID_EXPLICIT_RADIAL_FACTOR
