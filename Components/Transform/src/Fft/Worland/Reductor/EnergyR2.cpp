/**
 * @file EnergyR2.cpp
 * @brief Source of the implementation of the Worland R^2 energy operator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Reductor/EnergyR2.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"
#include "QuICC/SparseSM/Worland/R2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   EnergyR2::EnergyR2()
      : IWorlandEnergy(1)
   {
      this->mProfileId = Debug::Profiler::WORLANDREDU_ENERGYR2;
   }

   EnergyR2::~EnergyR2()
   {
   }

   void EnergyR2::makeOperator(Matrix& op, Matrix& eop, const internal::Array& igrid, const internal::Array& iweights, const Array& eweights, const int i) const
   {
      int l = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fastSize(i);

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;
      wnl.compute<MHDFloat>(op, nPoly, l, igrid, internal::Array(), ev::Set());

      nPoly = eweights.size();
      MHDFloat a = static_cast<MHDFloat>(wnl.alpha(2*l));
      MHDFloat b = static_cast<MHDFloat>(wnl.dBeta());
      ::QuICC::SparseSM::Worland::R2 spasm(nPoly, nPoly, a, b,  2*l);

      Matrix tmp(igrid.size(), nPoly);
      wnl.compute<MHDFloat>(tmp, nPoly, 2*l, igrid, iweights, ev::Set());
      eop.resize(igrid.size(), 1);
      eop = (eweights.transpose()*spasm.mat()*tmp.transpose()).transpose();
   }

   void EnergyR2::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      #if defined QUICC_WORLAND_REDUIMPL_MATRIX

         Matrix tmp = (this->mPOps.at(i)*in).array().abs2();
         rOut.transpose() = this->mEOps.at(i).transpose()*tmp;

      #elif defined QUICC_WORLAND_REDUIMPL_OTF
         int l = this->mspSetup->slow(i);
         int nPoly = this->mspSetup->fastSize(i);

         internal::Array igrid, iweights;
         this->computeEnergyQuadrature(igrid, iweights, this->mEWeights, i);
         this->mGrid = igrid.cast<MHDFloat>();
         this->mWeights = iweights.cast<MHDFloat>();

         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl bwnl;
         MatrixZ tmp(this->mGrid.size(), in.cols());
         bwnl.compute<MHDComplex>(tmp, nPoly, l, this->mGrid, internal::Array(), ev::OuterProduct<MHDComplex>(in));
         Matrix tmpC = tmp.array().abs2();

         nPoly = this->mEWeights.size();
         Polynomial::Worland::Wnl fwnl;
         Matrix tmpB(nPoly, in.cols());
         fwnl.compute<MHDFloat>(tmpB, nPoly, 2*l, this->mGrid, this->mWeights, ev::InnerProduct<MHDFloat>(tmpC));
         ::QuICC::SparseSM::Worland::R2 spasm(nPoly, nPoly, fwnl.alpha(2*l), fwnl.dBeta(), 2*l);
         rOut.transpose() = this->mEWeights.transpose()*spasm.mat()*tmpB;
      #endif //defined QUICC_WORLAND_REDUIMPL_MATRIX
   }

}
}
}
}
}
