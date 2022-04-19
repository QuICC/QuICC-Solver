/**
 * @file ITransformOperator.cpp
 * @brief Source of the interface for a generic transform operator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/ITransformOperator.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

   ITransformOperator::ITransformOperator()
      : mProfileId(Debug::Profiler::BLACKHOLE), mIsInitialized(false)
   {
   }

   ITransformOperator::~ITransformOperator()
   {
   }

   bool ITransformOperator::isInitialized() const
   {
      return this->mIsInitialized;
   }

   MHDFloat ITransformOperator::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}
