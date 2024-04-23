/**
 * @file IBesselRadialPower.cpp
 * @brief Source of the interface to a spherical Bessel radial grid power operator (e.g. energy)
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Reductor/IBesselRadialPower.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   IBesselRadialPower::IBesselRadialPower()
   {
   }

   void IBesselRadialPower::initOperators(const Internal::Array& igrid, const Internal::Array& iweights) const
   {
      // Reserve storage for the operators
      this->mOps.reserve(this->mspSetup->slowSize());

      // Loop over harmonic degrees
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         // Build operator
         this->mOps.push_back(Matrix(this->mspSetup->fastSize(i), igrid.size()));
         Matrix op;
         this->makeOperator(op, igrid, iweights, i);
         this->mOps.back() = op.transpose();
      }
   }

   void IBesselRadialPower::applyOperators(Matrix& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix(this->mProfileTag);

      // assert right sizes for input  matrix
      assert(in.cols() == this->mspSetup->blockSize());
      // assert right sizes for output matrix
      assert(rOut.rows() == this->outRows());
      assert(rOut.cols() == this->outCols());

      int start = 0;
      int outRows = this->mspSetup->fwdSize();
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         int cols = this->mspSetup->mult(i);
         int inRows = this->mspSetup->fastSize(i);

         this->applyOperator(rOut.block(0,start, outRows, cols), i, in.block(0,start, inRows, cols));

         start += cols;
      }
   }

   void IBesselRadialPower::applyOperators(MatrixZ& rOut, const MatrixZ& in) const
   {
      throw std::logic_error("Unused interface");
   }

   void IBesselRadialPower::defaultApplyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      rOut = (this->mOps.at(i).transpose()*in).array().abs2();
   }

   int IBesselRadialPower::outRows() const
   {
      return this->mspSetup->fwdSize();
   }

   int IBesselRadialPower::outCols() const
   {
      return this->mspSetup->blockSize();
   }

} // namespace Reductor
} // namespace Bessel
} // namespace Poly
} // Transform
} // QuICC
