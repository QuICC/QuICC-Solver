/**
 * @file IWorlandEnergy.cpp
 * @brief Source of the interface for a generic FFT based Worland energy reductor
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Reductor/IWorlandEnergy.hpp"

// Project includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Reduce.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   IWorlandEnergy::IWorlandEnergy(const int shift)
      : mcShift(shift)
   {
   }

   IWorlandEnergy::~IWorlandEnergy()
   {
   }

   void IWorlandEnergy::initBackend() const
   {
      #if defined QUICC_WORLAND_REDUIMPL_MATRIX
      // Reserve storage for the operators
      this->mPOps.reserve(this->mspSetup->slowSize());
      this->mEOps.reserve(this->mspSetup->slowSize());

      // Loop over harmonic degrees
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         // Energy calculation requires a different quadrature
         internal::Array igrid, iweights;
         Array eweights;
         this->computeEnergyQuadrature(igrid, iweights, eweights, i);

         // Build operator
         this->mPOps.push_back(Matrix(igrid.size(), this->mspSetup->fastSize(i)));
         this->mEOps.push_back(Matrix(igrid.size(), 1));
         this->makeOperator(this->mPOps.back(), this->mEOps.back(), igrid, iweights, eweights, i);
      }
      #elif defined QUICC_WORLAND_REDUIMPL_OTF

      internal::Array igrid, iweights;
      this->computeEnergyQuadrature(igrid, iweights, this->mEWeights, this->mspSetup->slowSize()-1);
      this->mGrid = igrid.cast<MHDFloat>();
      this->mWeights = iweights.cast<MHDFloat>();

      #endif //defined QUICC_WORLAND_REDUIMPL_MATRIX
   }

   void IWorlandEnergy::transform(Matrix& rOut, const MatrixZ& in) const
   {
      ProfilerMacro_start(Debug::Profiler::WORLANDTRA);
      ProfilerMacro_start(Debug::Profiler::WORLANDREDU);
      ProfilerMacro_start(this->mProfileId);

      assert(this->isInitialized());

      this->applyOperators(rOut, in);

      ProfilerMacro_stop(this->mProfileId);
      ProfilerMacro_stop(Debug::Profiler::WORLANDREDU);
      ProfilerMacro_stop(Debug::Profiler::WORLANDTRA);
   }

   void IWorlandEnergy::transform(Matrix& rOut, const Matrix& in) const
   {
      ProfilerMacro_start(Debug::Profiler::WORLANDTRA);
      ProfilerMacro_start(Debug::Profiler::WORLANDREDU);
      ProfilerMacro_start(this->mProfileId);

      assert(this->isInitialized());

      this->applyOperators(rOut, in);

      ProfilerMacro_stop(this->mProfileId);
      ProfilerMacro_stop(Debug::Profiler::WORLANDREDU);
      ProfilerMacro_stop(Debug::Profiler::WORLANDTRA);
   }

   void IWorlandEnergy::computeEnergyQuadrature(internal::Array& igrid, internal::Array& iweights, Array& eweights, const int iL) const
   {
      // Grid size for full integration of quadratic energy term (2N + L + 1) = 4/3 (3/2 N + 3/4L + 1)
      // and add shift
      int nrgSize = std::ceil(2.0*this->mspSetup->fastSize(iL) + this->mspSetup->slow(iL) + 1 + this->mcShift);

      Polynomial::Quadrature::WorlandRule quad;
      quad.computeQuadrature(igrid, iweights, nrgSize);

      // Create grids for energy calculations
      internal::Array ilegGrid, ilegWeights;
      Polynomial::Quadrature::LegendreRule lquad;
      lquad.computeQuadrature(ilegGrid, ilegWeights, nrgSize);
      ilegGrid.array() = ((ilegGrid.array() + 1.0)/2.0);

      int energySize = 2*this->mspSetup->fastSize(iL)+this->mcShift;
      int l = this->mspSetup->slow(iL);

      eweights.resize(energySize,1);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl wnl;
      wnl.compute<MHDFloat>(eweights, energySize, 2*l, ilegGrid, ilegWeights, ev::Reduce());
      eweights.array() *= 0.5;
   }

   void IWorlandEnergy::applyOperators(Matrix& rOut, const MatrixZ& in) const
   {
      // assert right sizes for input  matrix
      assert(in.cols() == this->mspSetup->blockSize());
      // assert right sizes for output matrix
      assert(rOut.rows() == this->mspSetup->blockSize());
      assert(rOut.cols() == 1);

      int start = 0;
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         int cols = this->mspSetup->mult(i);
         int inRows = this->mspSetup->specSize();

         this->applyOperator(rOut.block(start,0, cols, 1), i, in.block(0,start, inRows, cols));

         start += cols;
      }
   }

   int IWorlandEnergy::outRows() const
   {
      return this->mspSetup->blockSize();
   }

   int IWorlandEnergy::outCols() const
   {
      return 1;
   }

   MHDFloat IWorlandEnergy::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += IWorlandOperator::requiredStorage();

      // Storage for the operators
      for(auto it = this->mPOps.cbegin(); it != this->mPOps.cend(); ++it)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(it->size());
      }

      for(auto it = this->mEOps.cbegin(); it != this->mEOps.cend(); ++it)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(it->size());
      }

      // Storage for grid and weights
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(this->mGrid.size());
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(this->mWeights.size());
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(this->mEWeights.size());
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}
}
}
}
