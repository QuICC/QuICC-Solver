/**
 * @file I6DivR1D1R1_I4.cpp
 * @brief Source of the implementation of the Worland 1/R1 D R1 integrator but 0 mode is I4 P integrator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/I6DivR1D1R1_I4.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"
#include "QuICC/SparseSM/Worland/I6.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   void I6DivR1D1R1_I4<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      namespace ev = Polynomial::Worland::Evaluator;
      if(l == 0)
      {
         l = 1;
         Polynomial::Worland::Wnl wnl;

         // Internal computation uses dealiased modes
         const int extraN = 6*(!this->mcTruncQI); // I4 has 6 superdiagonals
         int nN = nPoly + extraN;
         this->checkGridSize(nN, l, igrid.size());

         Internal::Matrix tOp(igrid.size(), nN);

         wnl.compute<Internal::MHDFloat>(tOp, nN, l, igrid, iweights, ev::Set());

         MHDFloat a = static_cast<MHDFloat>(wnl.alpha(l));
         MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
         ::QuICC::SparseSM::Worland::I4 spasm(nN, nN, a, b, l);
         tOp = (spasm.mat()*tOp.transpose()).transpose();
         op = tOp.cast<MHDFloat>().leftCols(nPoly);
      }
      else
      {
         Polynomial::Worland::Wnl wnl;

         // Internal computation uses dealiased modes
         const int extraN = 9; // I6 has 9 superdiagonals
         int nN = nPoly + extraN;
         this->checkGridSize(nN, l, igrid.size());

         Internal::Matrix tOp(igrid.size(), nN);

         using poly_t = QuICC::Polynomial::Worland::implicit_t;
         Polynomial::Worland::r_1drWnl<poly_t> r_1drWnl;
         r_1drWnl.compute<Internal::MHDFloat>(tOp, nN, l, igrid, iweights, ev::Set());

         // Multiply by Quasi-inverse
         auto a = wnl.alpha(l);
         auto b = wnl.dBeta();
         ::QuICC::SparseSM::Worland::I6 spasm(nN, nN, a, b, l);
         tOp = (spasm.mat()*tOp.transpose()).transpose();
         op = tOp.cast<MHDFloat>().leftCols(nPoly);

         assert(op.rows() == igrid.size());
         assert(op.cols() == nPoly);
      }
   }

   void I6DivR1D1R1_I4<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
