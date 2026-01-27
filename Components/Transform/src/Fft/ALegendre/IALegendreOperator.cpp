/**
 * @file IALegendreOperator.cpp
 * @brief Source of the interface for a generic ALegendre FFT based operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/ALegendre/IALegendreOperator.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace ALegendre {

   IALegendreOperator::IALegendreOperator()
   {
      this->mProfileTag += "ALegendre::Fft";
   }

   IALegendreOperator::~IALegendreOperator()
   {
   }

   void IALegendreOperator::init(SharedTransformSetup spSetup) const
   {
      // Store the shared pointer to setup object
      this->mspSetup = std::dynamic_pointer_cast<IALegendreOperator::SetupType>(spSetup);

      //
      this->initBase();
   }

   void IALegendreOperator::init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const
   {
      throw std::logic_error("Unused interface");
   }

   void IALegendreOperator::transform(MatrixZ&, const MatrixZ&) const
   {
      throw std::logic_error("ALegendre FFT operator does not define a complex to complex transform");
   }

   void IALegendreOperator::transform(Matrix&, const MatrixZ&) const
   {
      throw std::logic_error("ALegendre FFT operator does not define a complex to real transform");
   }

   void IALegendreOperator::transform(Matrix&, const Matrix&) const
   {
      throw std::logic_error("ALegendre FFT operator does not define a real to real transform");
   }

   void IALegendreOperator::transform(MatrixZ&, const Matrix&) const
   {
      throw std::logic_error("ALegendre FFT operator does not define a real to complex transform");
   }

}
}
}
}
