/**
 * @file I4DivR1D1R1_Zero.cpp
 * @brief Source of the implementation of the Worland I4DivR1D1R1 integrator and l = 0 mode
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/BreakPoint.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/I4DivR1D1R1_Zero.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   I4DivR1D1R1_Zero::I4DivR1D1R1_Zero()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_I4DIVR1D1R1_ZERO;
   }

   I4DivR1D1R1_Zero::~I4DivR1D1R1_Zero()
   {
   }

   void I4DivR1D1R1_Zero::initBackend() const
   {
      // l = 0 mode is set to sero
      std::set<int> filter = {0};
      this->mBackend.setZFilter(filter);

      int lshift = 1; // operator shifts l by one
      int extraN = 6; // 6 extra modes are required due to I4 multiplication
      this->mBackend.init(*this->mspSetup, lshift, extraN);
      this->mBackend.addStorage(0, 1);
   }

   void I4DivR1D1R1_Zero::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      DivR1D1R1::computeWorlandExpansion(isEven);
      this->mBackend.applyI4(isEven);
      this->mBackend.output(rOut, isEven);
   }

   void I4DivR1D1R1_Zero::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      DivR1D1R1::computeWorlandExpansion(isEven);
      this->mBackend.applyI4(isEven);
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
