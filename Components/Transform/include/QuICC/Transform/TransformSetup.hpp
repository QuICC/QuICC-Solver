/** 
 * @file TransformSetup.hpp
 * @brief Implementation of base class for a generalized transform setup
 */

#ifndef QUICC_TRANSFORM_TRANSFORMSETUP_HPP
#define QUICC_TRANSFORM_TRANSFORMSETUP_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/GridPurpose.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of base class for a generalized transform setup
    */ 
   class TransformSetup
   {
      public:
         /**
          * @brief Constructor
          *
          * @param size       Size of the transform
          * @param specSize   Spectral output size (i.e without the padding)
          */
         TransformSetup(const int size, const int specSize, const GridPurpose::Id purpose);

         /**
          * @brief Constructor
          *
          * @param size       Size of the transform
          * @param blockSize  Number of similar transforms
          * @param specSize   Spectral output size (i.e without the padding)
          */
         TransformSetup(const int size, const int blockSize, const int specSize, const GridPurpose::Id purpose);

         /**
          * @brief Empty destructor
          */
         virtual ~TransformSetup();

         /**
          * @brief Purpose of grid values in transform
          */
         GridPurpose::Id purpose() const;

         /**
          * @brief Add index with multiplicity and array of fast indexes
          */
         void addIndex(const int slowIdx, const int mult, const ArrayI& fastIdx);

         /**
          * @brief Add index with multiplicity and size of fast indexes (0..size-1)
          */
         void addIndex(const int slowIdx, const int mult, const int size);

         /**
          * @brief Add index with multiplicity and assume full fast index list
          */
         virtual void addIndex(const int slowIdx, const int mult = 1);

         /**
          * @brief Get the size of the transform
          */
         int fwdSize() const;

         /**
          * @brief Get the total number of similar transforms
          */
         int blockSize() const;

         /**
          * @brief Get the spectral size of the transform
          */
         int specSize() const;

         /**
          * @brief Set box size scaling factor
          */
         virtual void setBoxScale(const MHDFloat boxScale);

         /**
          * @brief Lock setup to forbid adding new indexes
          */
         virtual void lock();

         /**
          * @brief Get the number of fast indexes
          */
         int fastSize(const int i) const;

         /**
          * @brief Get the number of slow indexes
          */
         int slowSize() const;

         /**
          * @brief Get the fast indexes
          */
         int fast(const int i, const int j) const;

         /**
          * @brief Get the slow indexes
          */
         int slow(const int i) const;

         /**
          * @brief Get the multipliers
          */
         int mult(const int i) const;
         
      protected:
         /**
          * @brief Purpose of grid values in transform
          */
         GridPurpose::Id mPurpose;

         /**
          * @brief Has full description for indexes
          */
         bool mHasFullIndexes;

         /**
          * @brief Size of the forward transform
          */
         int mFwdSize;

         /**
          * @brief Number of similar transforms
          */
         int mBlockSize;

         /**
          * @brief Spectral size
          */
         int mSpecSize;

         /**
          * @brief Is setup locked?
          */
         bool mIsLocked;
         /**
          * @brief Storage for the fast indexes
          */
         std::vector<ArrayI> mFast;

         /**
          * @brief Storage for the slow indexes
          */
         std::vector<int> mSlow;

         /**
          * @brief Storage for the multipliers
          */
         std::vector<int> mMult;

      private:

   };

   /// Typedef for an smart reference counting pointer for a TransformSetup
   typedef std::shared_ptr<TransformSetup>   SharedTransformSetup;

}
}

#endif // QUICC_TRANSFORM_TRANSFORMSETUP_HPP
