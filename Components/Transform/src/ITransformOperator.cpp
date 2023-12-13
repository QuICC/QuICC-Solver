/**
 * @file ITransformOperator.cpp
 * @brief Source of the interface for a generic transform operator
 */

// System includes
//
#include <cassert>
#include <boost/core/demangle.hpp>
#include <typeinfo>

// Project includes
//
#include "QuICC/Transform/ITransformOperator.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

   ITransformOperator::ITransformOperator()
      : mProfileTag(""), mIsInitialized(false)
   {
   }

   std::string ITransformOperator::opName() const
   {
      // \todo temporary fix, this needs cleaning
      std::string full = boost::core::demangle(typeid(*this).name());

      return full;
   }

   void ITransformOperator::setProfileTag()
   {
      this->mProfileTag = this->opName();
   }

   bool ITransformOperator::isInitialized() const
   {
      return this->mIsInitialized;
   }

   void ITransformOperator::init(SharedTransformSetup spSetup, const Internal::Array &igrid,
             const Internal::Array &iweights) const
   {
      throw std::logic_error("init needs to be implemented by the derived class");
   }

   void ITransformOperator::init(SharedTransformSetup spSetup) const
   {
      throw std::logic_error("init needs to be implemented by the derived class");
   }

   void ITransformOperator::transform(MatrixZ &rOut, const MatrixZ &in) const {}

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
