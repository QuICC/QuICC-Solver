/**
 * @file TransformSetup.cpp
 * @brief Source of basic transform setup class
 */

// System includes
//

// Class include
//
#include "QuICC/Transform/TransformSetup.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

   TransformSetup::TransformSetup(const int size, const int specSize, const GridPurpose::Id purpose)
      : mPurpose(purpose), mHasFullIndexes(true), mFwdSize(size), mBlockSize(-1), mSpecSize(specSize), mIsLocked(false)
   {
   }

   TransformSetup::TransformSetup(const int size, const int blockSize, const int specSize, const GridPurpose::Id purpose)
      : mPurpose(purpose), mHasFullIndexes(false), mFwdSize(size), mBlockSize(blockSize), mSpecSize(specSize), mIsLocked(false)
   {
   }

   TransformSetup::~TransformSetup()
   {
   }

   GridPurpose::Id TransformSetup::purpose() const
   {
      return this->mPurpose;
   }

   void TransformSetup::addIndex(const int slowIdx, const int mult, const ArrayI& fastIdx)
   {
      // Forbid modifying locked setup
      if(this->mIsLocked)
      {
         throw std::logic_error("Tried to add index to locked setup!");
      }

      this->mSlow.push_back(slowIdx);
      this->mMult.push_back(mult);
      this->mFast.push_back(fastIdx);
   }

   void TransformSetup::addIndex(const int slowIdx, const int mult, const int size)
   {
      ArrayI fastIdx = ArrayI::LinSpaced(size, 0, size-1);
      this->addIndex(slowIdx, mult, fastIdx);
   }

   void TransformSetup::addIndex(const int slowIdx, const int mult)
   {
      this->addIndex(slowIdx, mult, this->specSize());
   }

   int TransformSetup::fwdSize() const
   {
      return this->mFwdSize;
   }

   int TransformSetup::blockSize() const
   {
      return this->mBlockSize;
   }

   int TransformSetup::specSize() const
   {
      return this->mSpecSize;
   }

   void TransformSetup::setBoxScale(const MHDFloat)
   {
   }

   void TransformSetup::lock()
   {
      if(!this->mIsLocked)
      {
         // Only compute block size if not yet set
         if(this->mHasFullIndexes)
         {
            // Compute blockSize
            this->mBlockSize = 0;
            for(auto it = this->mMult.cbegin(); it != this->mMult.cend(); ++it)
            {
               this->mBlockSize += *it;
            }
         }

         this->mIsLocked = true;
      }
   }

   int TransformSetup::fastSize(const int i) const
   {
      return this->mFast.at(i).size();
   }

   int TransformSetup::slowSize() const
   {
      return this->mSlow.size();
   }

   int TransformSetup::fast(const int i, const int j) const
   {
      return this->mFast.at(j)(i);
   }

   int TransformSetup::slow(const int i) const
   {
      return this->mSlow.at(i);
   }

   int TransformSetup::mult(const int i) const
   {
      return this->mMult.at(i);
   }

}
}
