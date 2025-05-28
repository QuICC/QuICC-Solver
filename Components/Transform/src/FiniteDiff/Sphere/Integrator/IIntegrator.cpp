/**
 * @file IIntegrator.cpp
 * @brief Source of the interface to a Finite Difference sphere integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/FiniteDiff/Sphere/Integrator/IIntegrator.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace FiniteDiff {

namespace Sphere {

namespace Integrator {

   IIntegrator::IIntegrator() :
    IOperator()
   {
      this->mProfileTag += "-Integrator";
   }

   void IIntegrator::initOperators(const Internal::Array& igrid) const
   {
      // Reserve storage for the operators
      this->mOps.reserve(this->mspSetup->slowSize());

      // Loop over harmonic degrees
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         // Build operator
         this->mOps.push_back(Matrix(igrid.size(), this->mspSetup->fastSize(i)));
         this->makeOperator(this->mOps.back(), igrid, i);
      }
   }

   void IIntegrator::applyOperators(MatrixZ& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix(this->mProfileTag + "::applyOperators");

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

   void IIntegrator::applyOperators(Matrix& rOut, const MatrixZ& in) const
   {
      throw std::logic_error("Interface not used");
   }

   MHDFloat IIntegrator::requiredStorage() const
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

   void IIntegrator::defaultApplyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const
   {
      rOut = this->mOps.at(i).transpose()*in;
   }

   int IIntegrator::outRows() const
   {
      return this->mspSetup->fastSize(0);
   }

   int IIntegrator::outCols() const
   {
      return this->mspSetup->blockSize();
   }

} // namespace Integrator
} // namespace Sphere
} // namespace FiniteDiff
} // namespace Transform
} // namespace QuICC
