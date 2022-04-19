/**
 * @file Setup.cpp
 * @brief Source of Chebyshev FFT setup class
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

   Setup::Setup(const int size, const int blockSize, const int specSize, const GridPurpose::Id purpose)
      : ::QuICC::Transform::Fft::Setup(size, blockSize, specSize, purpose), mLower(std::numeric_limits<MHDFloat>::quiet_NaN()), mUpper(std::numeric_limits<MHDFloat>::quiet_NaN())
   {
      this->setBoxScale(1.0);
   }

   Setup::~Setup()
   {
   }

   void Setup::setBounds(const MHDFloat lower, const MHDFloat upper)
   {
      assert(lower < upper);

      this->mLower = lower;
      this->mUpper = upper;
   }

   void Setup::setBackwardSize()
   {
      this->mBwdSize = this->mFwdSize;
   }

   MHDFloat Setup::lower() const
   {
      assert(!std::isnan(this->mLower));

      return this->mLower;
   }

   MHDFloat Setup::upper() const
   {
      assert(!std::isnan(this->mLower));

      return this->mUpper;
   }

}
}
}
}
