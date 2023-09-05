/**
 * @file I4DivR1_Zero.cpp
 * @brief Source of the implementation of the Worland I4DivR1 integrator and l = 0 mode
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Integrator/Base/I4DivR1_Zero.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   I4DivR1_Zero<base_t>::I4DivR1_Zero()
   {
      this->setProfileTag();
   }


   void I4DivR1_Zero<base_t>::initBackend() const
   {
      // l = 0 mode is set to sero
      std::set<int> filter = {0};
      this->mBackend.setZFilter(filter);

      int lshift = -1; // operator shifts l by -1
      int extraN = 7; // 6 extra modes are required due to I4 multiplication, 1 for l - 1 shift
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

   void I4DivR1_Zero<base_t>::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      DivR1<base_t>::computeWorlandExpansion(isEven);
      this->mBackend.applyI4(isEven);
      this->mBackend.output(rOut, isEven);
   }

   void I4DivR1_Zero<base_t>::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      DivR1<base_t>::computeWorlandExpansion(isEven);
      this->mBackend.applyI4(isEven);
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
