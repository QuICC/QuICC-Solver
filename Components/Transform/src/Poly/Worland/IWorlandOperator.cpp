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
#include "QuICC/Transform/Poly/Worland/IWorlandOperator.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

   IWorlandOperator::IWorlandOperator()
      : ITransformOperator()
   {
      this->mProfileTag = "Worland::Poly";
   }

   IWorlandOperator::~IWorlandOperator()
   {
   }

   void IWorlandOperator::init(SharedTransformSetup spSetup, const internal::Array& igrid, const internal::Array& iweights) const
   {
      // Store the shared pointer to setup object
      if(spSetup)
      {
         this->mspSetup = std::dynamic_pointer_cast<IWorlandOperator::SetupType>(spSetup);
      } else
      {
         throw std::logic_error("Setup object is not initialized!");
      }

      // Initialise the operators
      this->initOperators(igrid, iweights);

      // Set initialization flag
      this->mIsInitialized = true;
   }

   void IWorlandOperator::init(SharedTransformSetup spSetup) const
   {
      throw std::logic_error("Unused interface");
   }

   void IWorlandOperator::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      assert(this->isInitialized());

      this->applyOperators(rOut, in);
   }

   void IWorlandOperator::transform(Matrix& rOut, const MatrixZ& in) const
   {
      assert(this->isInitialized());

      this->applyOperators(rOut, in);
   }

   void IWorlandOperator::applyOperators(MatrixZ&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Worland operator");
   }

   void IWorlandOperator::applyOperators(Matrix&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Worland operator");
   }

   void IWorlandOperator::checkGridSize(const int n, const int l, const int gN) const
   {
      int allowedN = (2*gN/3 - (l+1)/2 + 2);
      bool notValid = (n > allowedN);
      if(notValid)
      {
         throw std::logic_error("Worland grid is too small! (" + std::to_string(n) + " > " + std::to_string(allowedN) + ", n = " + std::to_string(n) + ", l = " + std::to_string(l) + ", gN = " + std::to_string(gN) + ")");
      }
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
