/**
 * @file IWorlandOperator.cpp
 * @brief Source of the interface for a Worland based transform operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Debug includes
//

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/BreakPoint.hpp"
#include "QuICC/Transform/Poly/Worland/IWorlandOperator.hpp"

// Project includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

   IWorlandOperator::IWorlandOperator()
      : ITransformOperator()
   {
   }

   IWorlandOperator::~IWorlandOperator()
   {
   }

   void IWorlandOperator::init(IWorlandOperator::SharedSetupType spSetup, const internal::Array& igrid, const internal::Array& iweights) const
   {
      // Store the shared pointer to setup object
      if(spSetup)
      {
         this->mspSetup = spSetup;
      } else
      {
         throw std::logic_error("Setup object is not initialized!");
      }

      // Initialise the operators
      this->initOperators(igrid, iweights);

      // Set initialization flag
      this->mIsInitialized = true;
   }

   void IWorlandOperator::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      ProfilerMacro_start(Debug::Profiler::WORLANDTRA);

      assert(this->isInitialized());

      this->applyOperators(rOut, in);

      ProfilerMacro_stop(Debug::Profiler::WORLANDTRA);
   }

   void IWorlandOperator::transform(Matrix& rOut, const MatrixZ& in) const
   {
      ProfilerMacro_start(Debug::Profiler::WORLANDTRA);

      assert(this->isInitialized());

      this->applyOperators(rOut, in);

      ProfilerMacro_stop(Debug::Profiler::WORLANDTRA);
   }

   void IWorlandOperator::applyOperators(MatrixZ&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Worland operator");
   }

   void IWorlandOperator::applyOperators(Matrix&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Worland operator");
   }

   MHDFloat IWorlandOperator::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += ITransformOperator::requiredStorage();
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}
}
}
