/**
 * @file IOperator.cpp
 * @brief Source of the interface for a Finite Differences based transform operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/FiniteDiff/Sphere/IOperator.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace FiniteDiff {

namespace Sphere {

   IOperator::IOperator()
      : ITransformOperator()
   {
      this->mProfileTag = "Sphere::FiniteDiff";
   }

   void IOperator::init(SharedTransformSetup spSetup, const Internal::Array& igrid) const
   {
      // Store the shared pointer to setup object
      if(spSetup)
      {
         this->mspSetup = std::dynamic_pointer_cast<IOperator::SetupType>(spSetup);
      } else
      {
         throw std::logic_error("Setup object is not initialized!");
      }

      // Initialise the operators
      this->initOperators(igrid);

      // Set initialization flag
      this->mIsInitialized = true;
   }

   void IOperator::init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const
   {
      throw std::logic_error("Unused interface");
   }

   void IOperator::init(SharedTransformSetup spSetup) const
   {
      throw std::logic_error("Unused interface");
   }

   void IOperator::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      assert(this->isInitialized());

      this->applyOperators(rOut, in);
   }

   void IOperator::transform(Matrix& rOut, const MatrixZ& in) const
   {
      assert(this->isInitialized());

      this->applyOperators(rOut, in);
   }

   void IOperator::applyOperators(MatrixZ&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Finite Differences operator");
   }

   void IOperator::applyOperators(Matrix&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Finite Differences operator");
   }

   void IOperator::checkGridSize(const int n, const int l, const int gN) const
   {
      int allowedN = (2*gN/3 - (l+1)/2 + 2);
      bool notValid = (n > allowedN);
      if(notValid)
      {
         throw std::logic_error("Finite Differences grid is too small! (" + std::to_string(n) + " > " + std::to_string(allowedN) + ", n = " + std::to_string(n) + ", l = " + std::to_string(l) + ", gN = " + std::to_string(gN) + ")");
      }
   }

   MHDFloat IOperator::requiredStorage() const
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
