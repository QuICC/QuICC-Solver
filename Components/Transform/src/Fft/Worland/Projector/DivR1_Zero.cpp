/**
 * @file DivR1_Zero.cpp
 * @brief Source of the implementation of the Worland DivR1_Zero projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Projector/DivR1_Zero.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   DivR1_Zero::DivR1_Zero()
   {
      this->mProfileId = Debug::Profiler::WORLANDPROJ_DIVR1_ZERO;
   }

   DivR1_Zero::~DivR1_Zero()
   {
   }

   void DivR1_Zero::initBackend() const
   {
      this->mBackend.init(*this->mspSetup, -1, true, true);
   }

   void DivR1_Zero::computeWorlandExpansion(const bool isEven) const
   {
      this->mBackend.scaleC(1.0/std::sqrt(Math::PI), isEven);
      this->mBackend.lowerBeta(-0.5, isEven);
      this->mBackend.backwardWorland(isEven);
   }

   void DivR1_Zero::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void DivR1_Zero::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->mBackend.output(rOut, isEven);
   }

   void DivR1_Zero::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void DivR1_Zero::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
