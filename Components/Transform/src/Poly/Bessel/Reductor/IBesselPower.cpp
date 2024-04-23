/**
 * @file IBesselPower.cpp
 * @brief Source of the interface to a spherical Bessel based reduction operator (e.g. energy)
 */

// System includes
//


// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Reductor/IBesselPower.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Quadrature/SphericalBesselRule.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   IBesselPower::IBesselPower(const int shift)
      : mcShift(shift)
   {
   }

   void IBesselPower::initOperators(const Internal::Array& icompgrid, const Internal::Array& icompweights) const
   {
      // Energy calculation requires a different quadrature
      Internal::Array igrid, iweights;
      this->computePowerQuadrature(igrid, iweights, icompgrid.size());

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
   }

   void IBesselPower::computePowerQuadrature(Internal::Array& igrid, Internal::Array& iweights, const int gSize) const
   {
      int nrgSize = gSize + 2*this->mcShift;

      Polynomial::Quadrature::SphericalBesselRule wquad;
      wquad.computeQuadrature(igrid, iweights, nrgSize);
   }

   void IBesselPower::applyOperators(Matrix& rOut, const MatrixZ& in) const
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

   void IBesselPower::applyOperators(MatrixZ& rOut, const MatrixZ& in) const
   {
      throw std::logic_error("Unused interface");
   }

   void IBesselPower::defaultApplyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      rOut = (this->mEOps.at(i).leftCols(rOut.rows()).transpose()*this->mOps.at(i)*in).array().abs2();
   }

   int IBesselPower::outRows() const
   {
      return this->mspSetup->fastSize(0);
   }

   int IBesselPower::outCols() const
   {
      return this->mspSetup->blockSize();
   }

}
}
}
}
}
