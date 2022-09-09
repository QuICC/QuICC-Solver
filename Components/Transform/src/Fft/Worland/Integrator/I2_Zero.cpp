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
      this->setProfileTag();
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

}
}
}
}
}
