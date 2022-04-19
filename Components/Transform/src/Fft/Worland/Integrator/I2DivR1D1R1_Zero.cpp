/**
 * @file I2DivR1D1R1_Zero.cpp
 * @brief Source of the implementation of the Worland I2DivR1D1R1 integrator and zero l = 0 mode
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/BreakPoint.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/I2DivR1D1R1_Zero.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   I2DivR1D1R1_Zero::I2DivR1D1R1_Zero()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_I2DIVR1D1R1_ZERO;
   }

   I2DivR1D1R1_Zero::~I2DivR1D1R1_Zero()
   {
   }

   void I2DivR1D1R1_Zero::initBackend() const
   {
      std::set<int> filter = {0};
      this->mBackend.setZFilter(filter);
      DivR1D1R1::initBackend();
   }

   void I2DivR1D1R1_Zero::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      DivR1D1R1::computeWorlandExpansion(isEven);
      this->mBackend.applyI2(isEven);
      this->mBackend.output(rOut, isEven);
   }

   void I2DivR1D1R1_Zero::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      DivR1D1R1::computeWorlandExpansion(isEven);
      this->mBackend.applyI2(isEven);
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
