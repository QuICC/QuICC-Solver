/**
 * @file P_Zero.cpp
 * @brief Source of the implementation of the Worland P_Zero integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Integrator/P_Zero.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   P_Zero::P_Zero()
      : IWorlandIntegrator()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_P;
   }

   P_Zero::~P_Zero()
   {
   }

   void P_Zero::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
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
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;

         // Internal computation uses dealiased modes
         int nN = nPoly;
         this->checkGridSize(nN, l, igrid.size());

         internal::Matrix tOp(igrid.size(), nN);

         wnl.compute<internal::MHDFloat>(tOp, nN, l, igrid, iweights, ev::Set());

         op = tOp.cast<MHDFloat>().leftCols(nPoly);
      }
   }

   void P_Zero::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_INTGIMPL_MATRIX
         this->defaultApplyOperator(rOut, i, in);
      #elif defined QUICC_WORLAND_INTGIMPL_OTF
         int l = this->mspSetup->slow(i);
         if(l == 0)
         {
            rOut.setZero();
         } else
         {
            int nPoly = this->mspSetup->fastSize(i);

            namespace ev = Polynomial::Worland::Evaluator;
            Polynomial::Worland::Wnl wnl;
            wnl.compute<MHDComplex>(rOut, nPoly, l, this->mGrid, this->mWeights, ev::InnerProduct<MHDComplex>(in));
         }
      #endif //defined QUICC_WORLAND_INTGIMPL_MATRIX
   }

}
}
}
}
}
