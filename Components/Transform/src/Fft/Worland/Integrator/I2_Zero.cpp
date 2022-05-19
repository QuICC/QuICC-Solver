/**
 * @file I2_Zero.cpp
 * @brief Source of the implementation of the Worland I2 integrator and zero l = 0 mode
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/BreakPoint.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/I2_Zero.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   I2_Zero::I2_Zero()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_I2_ZERO;
   }

   I2_Zero::~I2_Zero()
   {
   }

   void I2_Zero::initBackend() const
   {
      // l = 0 mode is set to sero
      std::set<int> filter = {0};
      this->mBackend.setZFilter(filter);

      int lshift = 0; // operator doesn't shift l
      int extraN = 3; // 3 extra modes are required due to I2 multiplication
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

   void I2_Zero::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      P::computeWorlandExpansion(isEven);
      this->mBackend.applyI2(isEven);
      this->mBackend.output(rOut, isEven);
   }

   void I2_Zero::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      P::computeWorlandExpansion(isEven);
      this->mBackend.applyI2(isEven);
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
