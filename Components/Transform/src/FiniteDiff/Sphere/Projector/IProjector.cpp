/**
 * @file IProjector.cpp
 * @brief Source of the interface to a Worland based projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/FiniteDiff/Sphere/Projector/IProjector.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace FiniteDiff {

namespace Sphere {

namespace Projector {

   IProjector::IProjector()
      : IOperator()
   {
      this->mProfileTag += "-Projector";
   }

   void IProjector::initOperators(const Internal::Array& igrid) const
   {
      // Reserve storage for the operators
      this->mOps.reserve(this->mspSetup->slowSize());

      // Loop over harmonic degrees
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         // Build operator
         this->mOps.push_back(Matrix(this->mspSetup->fastSize(i), igrid.size()));
         Matrix op;
         this->makeOperator(op, igrid, i);
         this->mOps.back() = op.transpose();
      }
   }

   void IProjector::applyOperators(MatrixZ& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix(this->mProfileTag + "::applyOperators");

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
         int inRows = this->mspSetup->fastSize(i);

         this->applyOperator(rOut.block(0, start, outRows, cols), i, in.block(0,start, inRows, cols));

         start += cols;
      }
   }

   void IProjector::applyOperators(Matrix& rOut, const MatrixZ& in) const
   {
      throw std::logic_error("Interface not used");
   }

   MHDFloat IProjector::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += IOperator::requiredStorage();

      // Storage for the operators
      for(auto it = this->mOps.cbegin(); it != this->mOps.cend(); ++it)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(it->size());
      }

      // Storage for grid and weights
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(this->mGrid.size());
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void IProjector::defaultApplyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const
   {
      rOut = this->mOps.at(i).transpose()*in;
   }

   int IProjector::outRows() const
   {
      return this->mspSetup->fwdSize();
   }

   int IProjector::outCols() const
   {
      return this->mspSetup->blockSize();
   }

}
}
}
}
}
