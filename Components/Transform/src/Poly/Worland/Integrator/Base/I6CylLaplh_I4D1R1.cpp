/**
 * @file I6CylLaplh_I4D1R1.cpp
 * @brief Source of the implementation of the Worland 1/R1 D R1 integrator but 0 mode is I4 D R1 integrator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/I6CylLaplh_I4D1R1.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/SparseSM/Worland/I4D1R1.hpp"
#include "QuICC/SparseSM/Worland/I6CylLaplh.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   void I6CylLaplh_I4D1R1<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;

      // Internal computation uses dealiased modes
      int nN = nPoly + 0;
      this->checkGridSize(nN, l, igrid.size());

      Internal::Matrix tOp(igrid.size(), nN);

      wnl.compute<Internal::MHDFloat>(tOp, nN, l, igrid, iweights, ev::Set());
      if(l == 0)
      {
         auto a = wnl.alpha(l);
         auto b = wnl.dBeta();
         ::QuICC::SparseSM::Worland::I4D1R1 spasm(nN, nN, a, b, 1);
         tOp = (spasm.mat()*tOp.transpose()).transpose();
      } else
      {
         auto a = wnl.alpha(l);
         auto b = wnl.dBeta();
         ::QuICC::SparseSM::Worland::I6CylLaplh spasm(nN, nN, a, b, l);
         tOp = (spasm.mat()*tOp.transpose()).transpose();
      }
      op = tOp.cast<MHDFloat>().leftCols(nPoly);
   }

   void I6CylLaplh_I4D1R1<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      this->defaultApplyOperator(rOut, i, in);
   }

}
}
}
}
}
