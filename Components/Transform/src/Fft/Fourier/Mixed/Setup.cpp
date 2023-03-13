/**
 * @file Setup.cpp
 * @brief Source of Mixed FFT setup class
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

   Setup::Setup(const int size, const int blockSize, const int specSize, const GridPurpose::Id purpose)
      : ::QuICC::Transform::Fft::Setup(size, blockSize, specSize, purpose)
   {
   }

   Setup::Setup(const int size, const int specSize, const GridPurpose::Id purpose)
      : ::QuICC::Transform::Fft::Setup(size, specSize, purpose)
   {
   }


   Setup::~Setup()
   {
   }

   void Setup::setBackwardSize()
   {
      this->mBwdSize = this->mFwdSize/2 + 1;
   }

}
}
}
}
}
