/**
 * @file I2.cpp
 * @brief Source of the implementation of the Worland I2 integrator and zero l = 0 mode
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Integrator/Base/I2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   I2<base_t>::I2()
   {
      this->setProfileTag();
   }


   void I2<base_t>::initBackend() const
   {
      int lshift = 0; // operator shifts l by one
      int extraN = 3; // 3 extra modes are required due to I2 multiplication
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

   void I2<base_t>::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      P<base_t>::computeWorlandExpansion(isEven);
      this->mBackend.applyI2(isEven);
      this->mBackend.output(rOut, isEven);
   }

   void I2<base_t>::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      P<base_t>::computeWorlandExpansion(isEven);
      this->mBackend.applyI2(isEven);
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
