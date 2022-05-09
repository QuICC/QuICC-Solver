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

   IALegendreOperator::~IALegendreOperator()
   {
   }

   void IALegendreOperator::init(IALegendreOperator::SharedSetupType spSetup, const internal::Array& igrid, const internal::Array& iweights) const
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Initialise the operators
      this->initOperators(igrid, iweights);

      // Set initialization flag
      this->mIsInitialized = true;
   }

   void IALegendreOperator::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix("IALegendreOperator::transformZ");
      assert(this->isInitialized());
      this->applyOperators(rOut, in);
   }

   void IALegendreOperator::transform(Matrix& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix("IALegendreOperator::transform");
      assert(this->isInitialized());
      this->applyOperators(rOut, in);
   }

   void IALegendreOperator::applyOperators(MatrixZ&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with ALegendre operator");
   }

   void IALegendreOperator::applyOperators(Matrix&, const MatrixZ&) const
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
