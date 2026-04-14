/**
 * @file IWorlandOperator.cpp
 * @brief Source of the interface for a generic Worland FFT based operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Worland/IWorlandOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

   IWorlandOperator::IWorlandOperator()
   {
      this->mProfileTag += "Worland::Fft";
   }

   void IWorlandOperator::init(SharedTransformSetup spSetup) const
   {
      // Store the shared pointer to setup object
      this->mspSetup = std::dynamic_pointer_cast<IWorlandOperator::SetupType>(spSetup);

      //
      this->initBase();
   }

   void IWorlandOperator::init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const
   {
      throw std::logic_error("Unused interface");
   }

   void IWorlandOperator::transform(Eigen::Ref<MatrixZ>, const Eigen::Ref<const MatrixZ>&) const
   {
      throw std::logic_error("Worland FFT operator does not define a complex to complex transform");
   }

   void IWorlandOperator::transform(Eigen::Ref<Matrix>, const Eigen::Ref<const MatrixZ>&) const
   {
      throw std::logic_error("Worland FFT operator does not define a complex to real transform");
   }

   void IWorlandOperator::transform(Eigen::Ref<Matrix>, const Eigen::Ref<const Matrix>&) const
   {
      throw std::logic_error("Worland FFT operator does not define a real to real transform");
   }

   void IWorlandOperator::transform(Eigen::Ref<MatrixZ>, const Eigen::Ref<const Matrix>&) const
   {
      throw std::logic_error("Worland FFT operator does not define a real to complex transform");
   }

} // namespace Worland
} // namespace Fft
} // namespace Transform
} // namespace QuICC
