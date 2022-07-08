/**
 * @file R1_Zero.cpp
 * @brief Source of the implementation of the Worland R1_Zero integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Integrator/R1_Zero.hpp"

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

   R1_Zero::R1_Zero()
      : IWorlandIntegrator()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_R1_ZERO;
   }

   R1_Zero::~R1_Zero()
   {
   }

   void R1_Zero::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
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
         Polynomial::Worland::Wnl wnl;

         // Internal computation uses dealiased modes
         int nN = nPoly;
         this->checkGridSize(nN, l, igrid.size());

         internal::Matrix tOp(igrid.size(), nN);

         wnl.compute<internal::MHDFloat>(tOp, nN, l, igrid, (iweights.array()*igrid.array()).matrix(), ev::Set());

         op = tOp.cast<MHDFloat>().leftCols(nPoly);
      }
   }

   void R1_Zero::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
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
            namespace ev = Polynomial::Worland::Evaluator;
            Polynomial::Worland::Wnl wnl;

            wnl.compute<MHDComplex>(rOut, this->mspSetup->fastSize(i), this->mspSetup->slow(i), this->mGrid, (this->mWeights.array()*this->mGrid.array()).matrix(), ev::InnerProduct<MHDComplex>(in));
         }
      #endif //defined QUICC_WORLAND_INTGIMPL_MATRIX
   }

}
}
}
}
}
