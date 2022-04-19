/**
 * @file R1.cpp
 * @brief Source of the implementation of the Worland R1 integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Integrator/R1.hpp"

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

   R1::R1()
      : IWorlandIntegrator()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_R1;
   }

   R1::~R1()
   {
   }

   void R1::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
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

         wnl.compute<MHDFloat>(op, nPoly, l, igrid, (iweights.array()*igrid.array()).matrix(), ev::Set());
      }
   }

   void R1::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
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
