/**
 * @file I6DivR1_Zero.cpp
 * @brief Source of the implementation of the Worland 1/R1 integrator but 0 mode is zeroed
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/I6DivR1_Zero.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"
#include "QuICC/SparseSM/Worland/I6.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   void I6DivR1_Zero<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      if(l == 0)
      {
         op.setZero();
      }
      else
      {
         // Internal computation uses dealiased modes
         const int extraN = 9; // I6 has 9 superdiagonals
         int nN = nPoly + extraN;
         this->checkGridSize(nN, l, igrid.size());

         Internal::Matrix tOp(igrid.size(), nN);

         using poly_t = QuICC::Polynomial::Worland::implicit_t;
         Polynomial::Worland::r_1Wnl<poly_t> r_1Wnl;
         r_1Wnl.compute<Internal::MHDFloat>(tOp, nN, l, igrid, iweights, ev::Set());

         // Multiply by Quasi-inverse
         auto a = r_1Wnl.alpha(l);
         auto b = r_1Wnl.dBeta();
         ::QuICC::SparseSM::Worland::I6 spasm(nN, nN, a, b, l);
         tOp = (spasm.mat()*tOp.transpose()).transpose();
         op = tOp.cast<MHDFloat>().leftCols(nPoly);

         assert(op.rows() == igrid.size());
         assert(op.cols() == nPoly);
      }
   }

   void I6DivR1_Zero<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
