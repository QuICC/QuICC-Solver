/**
 * @file I2.cpp
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
#include "QuICC/Transform/Fft/Worland/Integrator/I2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   I2::I2()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_I2;
   }

   I2::~I2()
   {
   }

   void I2::initBackend() const
   {
      int lshift = 0; // operator shifts l by one
      int extraN = 3; // 3 extra modes are required due to I2 multiplication
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

   void I2::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      P::computeWorlandExpansion(isEven);
      this->mBackend.applyI2(isEven);
      this->mBackend.output(rOut, isEven);
   }

   void I2::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
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
