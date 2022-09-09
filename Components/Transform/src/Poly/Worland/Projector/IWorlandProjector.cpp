/**
 * @file IWorlandProjector.cpp
 * @brief Source of the interface to a Worland based projector
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Projector/IWorlandProjector.hpp"

// Project includes
//
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   IWorlandProjector::IWorlandProjector()
      : IWorlandOperator()
   {
      this->mProfileTag += "-Projector";
   }

   IWorlandProjector::~IWorlandProjector()
   {
   }

   void IWorlandProjector::initOperators(const internal::Array& igrid, const internal::Array& iweights) const
   {
      #if defined QUICC_WORLAND_PROJIMPL_MATRIX
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

      #elif defined QUICC_WORLAND_PROJIMPL_OTF

         // Store grid and weights
         this->mGrid = igrid;
         this->mWeights = iweights;

      #endif //defined QUICC_WORLAND_PROJIMPL_MATRIX
   }

   void IWorlandProjector::applyOperators(MatrixZ& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix(this->mProfileTag);

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

   MHDFloat IWorlandProjector::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += IWorlandOperator::requiredStorage();

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

   void IWorlandProjector::defaultApplyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      rOut = this->mOps.at(i).transpose()*in;
   }

   int IWorlandProjector::outRows() const
   {
      return this->mspSetup->fwdSize();
   }

   int IWorlandProjector::outCols() const
   {
      return this->mspSetup->blockSize();
   }

}
}
}
}
}
