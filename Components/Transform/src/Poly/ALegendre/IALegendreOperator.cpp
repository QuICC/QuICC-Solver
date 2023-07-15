/**
 * @file IALegendreOperator.cpp
 * @brief Source of the interface for a associated Legendre based transform operator
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/IALegendreOperator.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

   IALegendreOperator::IALegendreOperator()
      : ITransformOperator()
   {
   }

   void IALegendreOperator::init(SharedTransformSetup spSetup, const OpArray& igrid, const OpArray& iweights) const
   {
      // Store the shared pointer to setup object
      this->mspSetup = std::dynamic_pointer_cast<IALegendreOperator::SetupType>(spSetup);

      // Initialise the operators
      this->initOperators(igrid, iweights);

      // Set initialization flag
      this->mIsInitialized = true;
   }

   void IALegendreOperator::init(SharedTransformSetup spSetup) const
   {
      throw std::logic_error("Unused interface");
   }

   void IALegendreOperator::transform(OpMatrixZ& rOut, const OpMatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix("IALegendreOperator::transformZ");
      assert(this->isInitialized());
      this->applyOperators(rOut, in);
   }

   void IALegendreOperator::transform(OpMatrix& rOut, const OpMatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix("IALegendreOperator::transform");
      assert(this->isInitialized());
      throw std::logic_error("Data is not compatible with ALegendre operator");
   }

   void IALegendreOperator::applyOperators(OpMatrixZ&, const OpMatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with ALegendre operator");
   }

   MHDFloat IALegendreOperator::requiredStorage() const
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
