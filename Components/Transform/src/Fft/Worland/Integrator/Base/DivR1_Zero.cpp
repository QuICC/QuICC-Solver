/**
 * @file DivR1_Zero.cpp
 * @brief Source of the implementation of the Worland DivR1_Zero integrator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Integrator/Base/DivR1_Zero.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   DivR1_Zero<base_t>::DivR1_Zero()
   {
      this->setProfileTag();
   }


   void DivR1_Zero<base_t>::initBackend() const
   {
      // l = 0 mode is set to sero
      std::set<int> filter = {0};
      this->mBackend.setZFilter(filter);

      int lshift = -1; // operator shifts l by one
      int extraN = 0; // no extra modes are required
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

}
}
}
}
}
