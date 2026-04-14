/**
 * @file Setup.cpp
 * @brief Source of Complex FFT setup class
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

   Setup::Setup(const int size, const int specSize, const GridPurpose::Id purpose)
      : ::QuICC::Transform::Fft::Setup(size, specSize, purpose)
   {
   }

   Setup::Setup(const int size, const int blockSize, const int specSize, const GridPurpose::Id purpose)
      : ::QuICC::Transform::Fft::Setup(size, blockSize, specSize, purpose)
   {
   }

   void Setup::lock()
   {
      if(!this->mIsLocked)
      {
         this->mIdBlocks.resize(this->slowSize(),2);

         for(int i = 0; i < this->slowSize(); ++i)
         {
            this->mIdBlocks(i,0) = this->slow(i);
            this->mIdBlocks(i,1) = this->mult(i);
         }

         ::QuICC::Transform::Fft::Setup::lock();
      }
   }

   void Setup::setBackwardSize()
   {
      this->mBwdSize = this->mFwdSize;
   }

   const MatrixI& Setup::idBlocks() const
   {
      assert(this->mIsLocked);

      return this->mIdBlocks;
   }

}
}
}
}
}
