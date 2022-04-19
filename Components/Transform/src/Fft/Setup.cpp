/**
 * @file Setup.cpp
 * @brief Source of FFT setup class
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Setup.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

   Setup::Setup(const int size, const int specSize, const GridPurpose::Id purpose)
      : TransformSetup(size, specSize, purpose), mBwdSize(0), mBoxScale(-4242)
   {
   }

   Setup::Setup(const int size, const int blockSize, const int specSize, const GridPurpose::Id purpose)
      : TransformSetup(size, blockSize, specSize, purpose), mBwdSize(0), mBoxScale(-4242)
   {
   }

   Setup::~Setup()
   {
   }

   void Setup::lock()
   {
      if(!this->mIsLocked)
      {
         this->setBackwardSize();

         // Safety assert
         assert(this->mBwdSize >= this->mSpecSize);
         assert(this->mFwdSize >= this->mBwdSize);
         assert(this->mBoxScale != -4242);

         TransformSetup::lock();
      }
   }

   void Setup::setBoxScale(const MHDFloat boxScale)
   {
      this->mBoxScale = boxScale;
   }

   int Setup::bwdSize() const
   {
      assert(this->mIsLocked);

      return this->mBwdSize;
   }

   int Setup::padSize() const
   {
      assert(this->mIsLocked);

      return this->mBwdSize - this->mSpecSize;
   }

   MHDFloat Setup::boxScale() const
   {
      return this->mBoxScale;
   }

}
}
}
