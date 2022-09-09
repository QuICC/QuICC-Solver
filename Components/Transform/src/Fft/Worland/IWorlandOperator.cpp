/**
 * @file IWorlandOperator.cpp
 * @brief Source of the interface for a generic Worland FFT based operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/IWorlandOperator.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

   IWorlandOperator::IWorlandOperator()
   {
      this->mProfileTag += "Worland::Fft";
   }

   IWorlandOperator::~IWorlandOperator()
   {
   }

   void IWorlandOperator::init(IWorlandOperator::SharedSetupType spSetup) const
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      //
      this->initBase();
   }

   void IWorlandOperator::transform(MatrixZ&, const MatrixZ&) const
   {
      throw std::logic_error("Worland FFT operator does not define a complex to complex transform");
   }

   void IWorlandOperator::transform(Matrix&, const MatrixZ&) const
   {
      throw std::logic_error("Worland FFT operator does not define a complex to real transform");
   }

   void IWorlandOperator::transform(Matrix&, const Matrix&) const
   {
      throw std::logic_error("Worland FFT operator does not define a real to real transform");
   }

   void IWorlandOperator::transform(MatrixZ&, const Matrix&) const
   {
      throw std::logic_error("Worland FFT operator does not define a real to complex transform");
   }

}
}
}
}
