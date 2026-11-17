/**
 * @file Setup.cpp
 * @brief Source of Worland FFT setup class
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

   Setup::Setup(const int size, const int specSize, const GridPurpose::Id purpose)
      : ::QuICC::Transform::Fft::Setup(size, specSize, purpose)
   {
      this->setBoxScale(1.0);
   }

   void Setup::setBackwardSize()
   {
      this->mBwdSize = this->mFwdSize;
   }

}
}
}
}
