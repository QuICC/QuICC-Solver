/**
 * @file IWorlandPower.cpp
 * @brief Source of the interface to a Worland based reduction operator (e.g. energy)
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandPower.hpp"

// Project includes
//
#include "Profiler/Interface.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandSphEnergyRule.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Reduce.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   IWorlandPower::IWorlandPower(const int shift)
      : mcShift(shift)
   {
   }

   IWorlandPower::~IWorlandPower()
   {
   }

   void IWorlandPower::initOperators(const internal::Array& icompgrid, const internal::Array& icompweights) const
   {
      // Energy calculation requires a different quadrature
      internal::Array igrid, iweights;
      this->computePowerQuadrature(igrid, iweights, icompgrid.size());

      #if defined QUICC_WORLAND_REDUIMPL_MATRIX
      // Reserve storage for the operators
      this->mOps.reserve(this->mspSetup->slowSize());
      this->mEOps.reserve(this->mspSetup->slowSize());

      // Loop over harmonic degrees
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         // Build operator
         this->mOps.push_back(Matrix(igrid.size(), this->mspSetup->fastSize(i)));
         this->mEOps.push_back(Matrix(igrid.size(), 1));
         this->makeOperator(this->mOps.back(), this->mEOps.back(), igrid, iweights, i);
      }
      #elif defined QUICC_WORLAND_REDUIMPL_OTF

      this->mGrid = igrid.cast<MHDFloat>();
      this->mWeights = iweights.cast<MHDFloat>();

      #endif //defined QUICC_WORLAND_REDUIMPL_MATRIX
   }

   void IWorlandPower::computePowerQuadrature(internal::Array& igrid, internal::Array& iweights, const int gSize) const
   {
      int nrgSize = gSize + 2*this->mcShift;

      Polynomial::Quadrature::WorlandSphEnergyRule wquad;
      wquad.computeQuadrature(igrid, iweights, nrgSize);
   }

   void IWorlandPower::applyOperators(Matrix& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix(this->mProfileTag);

      // assert right sizes for input  matrix
      assert(in.cols() == this->mspSetup->blockSize());
      // assert right sizes for output matrix
      assert(rOut.rows() == this->outRows());
      assert(rOut.cols() == this->outCols());

      int start = 0;
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         int cols = this->mspSetup->mult(i);
         int inRows = this->mspSetup->fastSize(i);

         this->applyOperator(rOut.block(0,start, inRows, cols), i, in.block(0,start, inRows, cols));

         start += cols;
      }
   }

   void IWorlandPower::defaultApplyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      rOut = (this->mEOps.at(i).leftCols(rOut.rows()).transpose()*this->mOps.at(i)*in).array().abs2();
   }

   int IWorlandPower::outRows() const
   {
      return this->mspSetup->fastSize(0);
   }

   int IWorlandPower::outCols() const
   {
      return this->mspSetup->blockSize();
   }

}
}
}
}
}
