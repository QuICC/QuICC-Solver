/**
 * @file EnergyD1R1.cpp
 * @brief Source of the implementation of the Worland D R energy operator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/BreakPoint.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/EnergyD1R1.hpp"

// Project includes
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/drWnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   EnergyD1R1::EnergyD1R1()
      : IWorlandEnergy(0)
   {
      this->mProfileId = Debug::Profiler::WORLANDREDU_ENERGYD1R1;
   }

   EnergyD1R1::~EnergyD1R1()
   {
   }

   void EnergyD1R1::makeOperator(Matrix& op, Matrix& eop, const internal::Array& igrid, const internal::Array& iweights, const Array& eweights, const int i) const
   {
      int l = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fastSize(i);

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::drWnl drwnl;
      drwnl.compute<MHDFloat>(op, nPoly, l, igrid, internal::Array(), ev::Set());

      nPoly = eweights.size();
      Matrix tmp(igrid.size(), nPoly);
      Polynomial::Worland::Wnl wnl;
      wnl.compute<MHDFloat>(tmp, nPoly, 2*l, igrid, iweights, ev::Set());
      eop.resize(igrid.size(), 1);
      eop = (eweights.transpose()*tmp.transpose()).transpose();
   }

   void EnergyD1R1::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
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
         Polynomial::Worland::drWnl bwnl;
         MatrixZ tmp(this->mGrid.size(), in.cols());
         bwnl.compute<MHDComplex>(tmp, nPoly, l, this->mGrid, internal::Array(), ev::OuterProduct<MHDComplex>(in));
         Matrix tmpC = tmp.array().abs2();

         nPoly = this->mEWeights.size();
         Polynomial::Worland::Wnl fwnl;
         Matrix tmpB(nPoly, in.cols());
         fwnl.compute<MHDFloat>(tmpB, nPoly, 2*l, this->mGrid, this->mWeights, ev::InnerProduct<MHDFloat>(tmpC));
         rOut.transpose() = this->mEWeights.transpose()*tmpB;
      #endif //defined QUICC_WORLAND_REDUIMPL_MATRIX
   }

}
}
}
}
}
