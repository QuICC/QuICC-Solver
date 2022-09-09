/**
 * @file ITransformOperator.cpp
 * @brief Source of the interface for a generic transform operator
 */

// System includes
//
#include <cassert>
#include <boost/core/demangle.hpp>
#include <typeinfo>

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
      : mProfileTag(""), mIsInitialized(false)
   {
   }

   ITransformOperator::~ITransformOperator()
   {
   }

   std::string ITransformOperator::opName() const
   {
      std::string full = boost::core::demangle(typeid(*this).name());
      std::size_t pos = full.rfind(':');
      std::string op = full.substr(pos+1, full.size()-pos);

      return op;
   }

   void ITransformOperator::setProfileTag()
   {
      this->mProfileTag += "-" + this->opName();
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
