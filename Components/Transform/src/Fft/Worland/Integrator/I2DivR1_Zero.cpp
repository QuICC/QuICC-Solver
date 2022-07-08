/**
 * @file I2DivR1_Zero.cpp
 * @brief Source of the implementation of the Worland I2DivR1 integrator and l = 0 zero
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/BreakPoint.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/I2DivR1_Zero.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   I2DivR1_Zero::I2DivR1_Zero()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_I2DIVR1_ZERO;
   }

   I2DivR1_Zero::~I2DivR1_Zero()
   {
   }

   void I2DivR1_Zero::initBackend() const
   {
      // l = 0 mode is set to sero
      std::set<int> filter = {0};
      this->mBackend.setZFilter(filter);

      int lshift = -1; // operator shifts l by -1
      int extraN = 4; // 3 extra modes are required due to I2 multiplication, 1 for l - 1 shift
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

   void I2DivR1_Zero::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      DivR1::computeWorlandExpansion(isEven);
      this->mBackend.applyI2(isEven);
      this->mBackend.output(rOut, isEven);
   }

   void I2DivR1_Zero::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      DivR1::computeWorlandExpansion(isEven);
      this->mBackend.applyI2(isEven);
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
