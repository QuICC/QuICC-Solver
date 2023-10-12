/**
 * @file IWorlandIntegrator.cpp
 * @brief Source of the interface to a Worland based integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Integrator/IWorlandIntegrator.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   IWorlandIntegrator::IWorlandIntegrator()
      : IWorlandOperator(),
#ifdef QUICC_TRANSFORM_WORLAND_TRUNCATE_QI
        mcTruncQI(true)
#else
        mcTruncQI(false)
#endif // QUICC_TRANSFORM_WORLAND_TRUNCATE_QI
   {
      this->mProfileTag += "-Integrator";
   }

   IWorlandIntegrator::~IWorlandIntegrator()
   {
   }

   void IWorlandIntegrator::initOperators(const Internal::Array& igrid, const Internal::Array& iweights) const
   {
      // Reserve storage for the operators
      this->mOps.reserve(this->mspSetup->slowSize());

      // Loop over harmonic degrees
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         // Build operator
         this->mOps.push_back(Matrix(igrid.size(), this->mspSetup->fastSize(i)));
         this->makeOperator(this->mOps.back(), igrid, iweights, i);
      }
   }

   void IWorlandIntegrator::applyOperators(MatrixZ& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix(this->mProfileTag);

      // assert right sizes for input matrix
      assert(in.rows() == this->mspSetup->fwdSize());
      assert(in.cols() == this->mspSetup->blockSize());
      // assert right sizes for output matrix
      assert(rOut.cols() == this->mspSetup->blockSize());

      int start = 0;
      int inRows = this->mspSetup->fwdSize();
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         int cols = this->mspSetup->mult(i);
         int outRows = this->mspSetup->fastSize(i);

         this->applyOperator(rOut.block(0, start, outRows, cols), i, in.block(0,start, inRows, cols));

         start += cols;
      }
   }

   void IWorlandIntegrator::applyOperators(Matrix& rOut, const MatrixZ& in) const
   {
      throw std::logic_error("Interface not used");
   }

   MHDFloat IWorlandIntegrator::requiredStorage() const
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

   void IWorlandIntegrator::defaultApplyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      rOut = this->mOps.at(i).transpose()*in;
   }

   int IWorlandIntegrator::outRows() const
   {
      return this->mspSetup->fastSize(0);
   }

   int IWorlandIntegrator::outCols() const
   {
      return this->mspSetup->blockSize();
   }

}
}
}
}
}
