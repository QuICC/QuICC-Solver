/**
 * @file DivR1D1R1_Zero.cpp
 * @brief Source of the implementation of the Worland DivR1D1R1_Zero integrator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Integrator/Base/DivR1D1R1_Zero.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   DivR1D1R1_Zero<base_t>::DivR1D1R1_Zero()
   {
      this->setProfileTag();
   }


   void DivR1D1R1_Zero<base_t>::initBackend() const
   {
      // l = 0 mode is set to sero
      std::set<int> filter = {0};
      this->mBackend.setZFilter(filter);

      int lshift = -1; // operator shifts l by one
      int extraN = 1; // 1 extra mode is required for l-1 shift
      this->mBackend.init(*this->mspSetup, lshift, extraN);
      this->mBackend.addStorage(0, 1);
   }

}
}
}
}
}
