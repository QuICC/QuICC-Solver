/**
 * @file IALegendreProjector.cpp
 * @brief Source of the interface to a associated Legendre based projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Projector/IALegendreProjector.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   IALegendreProjector::IALegendreProjector()
      : IALegendreOperator()
   {
   }

   void IALegendreProjector::initOperators(const OpArray& igrid, const OpArray& iweights) const
   {
      // Initit specialized data for operators
      this->initSpecial();

      // reserve storage
      this->mOps.reserve(this->mspSetup->slowSize());

      // Loop over harmonic orders
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         // Build operator
         this->mOps.push_back(Matrix(this->mspSetup->fastSize(i), igrid.size()));
         Matrix op;
         this->makeOperator(op, igrid, iweights, i);
         this->mOps.back() = op.transpose();
      }
   }

   void IALegendreProjector::applyOperators(OpMatrixZ& rOut, const OpMatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix("IALegendreProjector::applyOperators");

      // assert right sizes for input  matrix
      assert(in.cols() == this->mspSetup->blockSize());
      // assert right sizes for output matrix
      assert(rOut.rows() == this->mspSetup->fwdSize());
      assert(rOut.cols() == this->mspSetup->blockSize());

      int start = 0;
      int outRows = this->mspSetup->fwdSize();
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         int cols = this->mspSetup->mult(i);
         int m = this->mspSetup->slow(i);
         int inRows = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1;

         this->applyOperator(rOut.block(0, start, outRows, cols), i, in.block(0, start, inRows, cols));

         start += cols;
      }
   }

   void IALegendreProjector::initSpecial() const
   {
   }

   int IALegendreProjector::outRows() const
   {
      return this->mspSetup->fwdSize();
   }

   int IALegendreProjector::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   MHDFloat IALegendreProjector::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += IALegendreOperator::requiredStorage();

      // Storage for the operators
      for(auto it = this->mOps.cbegin(); it != this->mOps.cend(); ++it)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(it->size());
      }

      // Storage for grid and weights
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(this->mGrid.size());
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(this->mWeights.size());
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}
}
}
}
