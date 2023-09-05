/**
 * @file DivR1_Zero.cpp
 * @brief Source of the implementation of the Worland DivR1_Zero projector
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Projector/Base/DivR1_Zero.hpp"
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   DivR1_Zero<base_t>::DivR1_Zero()
   {
      this->setProfileTag();
   }

   void DivR1_Zero<base_t>::initBackend() const
   {
      int lshift = -1; // operator shifts l by -1
      int extraN = 0; // no modes are required
      bool onlyShiftParity = true;
      bool zeroNegativeL = true;
      this->mBackend.init(*this->mspSetup, lshift, extraN, onlyShiftParity, zeroNegativeL);
   }

   void DivR1_Zero<base_t>::computeWorlandExpansion(const bool isEven) const
   {
      this->mBackend.scaleC(1.0/std::sqrt(Math::PI), isEven);
      this->mBackend.lowerBeta(-0.5, isEven);
      this->mBackend.backwardWorland(isEven);
   }

   void DivR1_Zero<base_t>::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void DivR1_Zero<base_t>::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->mBackend.output(rOut, isEven);
   }

   void DivR1_Zero<base_t>::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void DivR1_Zero<base_t>::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
