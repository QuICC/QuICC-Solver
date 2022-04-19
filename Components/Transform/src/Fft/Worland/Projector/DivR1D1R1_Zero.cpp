/**
 * @file DivR1D1R1_Zero.cpp
 * @brief Source of the implementation of the Worland DivR1D1R1_Zero projector
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/BreakPoint.hpp"
#include "QuICC/Transform/Fft/Worland/Projector/DivR1D1R1_Zero.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   DivR1D1R1_Zero::DivR1D1R1_Zero()
   {
      this->mProfileId = Debug::Profiler::WORLANDPROJ_DIVR1D1R1_ZERO;
   }

   DivR1D1R1_Zero::~DivR1D1R1_Zero()
   {
   }

   void DivR1D1R1_Zero::initBackend() const
   {
      this->mBackend.init(*this->mspSetup, -1, true, true);
      this->mBackend.addStorage(1, 0);
      this->mBackend.setWSize(1);
   }

   void DivR1D1R1_Zero::computeWorlandExpansion(const bool isEven) const
   {
      const int main = 0;
      const int extra = 1;

      this->mBackend.scaleC(1.0/std::sqrt(Math::PI), isEven);

      // Copy expansion shifting indexes by -1 and scaling
      this->mBackend.copy(extra, main, -1, isEven);
      this->mBackend.scaleD(isEven, 0, extra);
      this->mBackend.lshift(extra, 1, isEven);
      this->mBackend.lowerAlpha(0.5, isEven, extra, 1.0);
      this->mBackend.lowerR2Beta(-0.5, isEven, extra, 1.0);

      // Scaling original expansion
      this->mBackend.scaleALPY(1.0, 1.0, isEven);

      // Add second term
      this->mBackend.add(main, extra, 0, isEven);
      this->mBackend.lowerBeta(-0.5, isEven);

      this->mBackend.backwardWorland(isEven);
   }

   void DivR1D1R1_Zero::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void DivR1D1R1_Zero::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->mBackend.output(rOut, isEven);
   }

   void DivR1D1R1_Zero::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void DivR1D1R1_Zero::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
