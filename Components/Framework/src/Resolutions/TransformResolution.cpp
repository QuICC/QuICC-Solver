/** 
 * @file TransformResolution.cpp
 * @brief Source of the resolution object for a transform
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Resolutions/TransformResolution.hpp"

// Project includes
//

namespace QuICC {

   TransformResolution::TransformResolution(const std::vector<ArrayI>& fwd, const std::vector<ArrayI>& bwd, const std::vector<ArrayI>& idx2D, const ArrayI& idx3D)
      : mFwd(fwd), mBwd(bwd), mIdx2D(idx2D), mIdx3D(idx3D), mDimF1D(idx3D.size()), mDimB1D(idx3D.size()), mDim2D(idx3D.size()), mDim3D(idx3D.size())
   {
      // Initialise the dimensions
      this->initDimensions();
   }

   TransformResolution::~TransformResolution()
   {
   }

   bool TransformResolution::isCleared() const
   {
      return (this->mIdx3D.size() == 0);
   }

   void TransformResolution::initDimensions()
   {
      for(int k = 0; k < this->mDim3D; k++)
      {
         this->mDimF1D(k) = this->mFwd.at(k).rows();
         this->mDimB1D(k) = this->mBwd.at(k).rows();
         this->mDim2D(k) = this->mIdx2D.at(k).size();
      }
   }

   template <> int TransformResolution::dim<Dimensions::Data::DATF1D>(const int k) const
   {
      assert(this->mDimF1D.size() > k);

      return this->mDimF1D(k);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DATF1D>() const
   {
      return this->dim<Dimensions::Data::DATF1D>(0);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DATB1D>(const int k) const
   {
      assert(this->mDimB1D.size() > k);

      return this->mDimB1D(k);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DATB1D>() const
   {
      return this->dim<Dimensions::Data::DATB1D>(0);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DAT2D>(const int k) const
   {
      assert(this->mDim2D.size() > k);

      return this->mDim2D(k);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DAT2D>() const
   {
      return this->dim<Dimensions::Data::DAT2D>(0);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DAT3D>() const
   {
      return this->mDim3D;
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATF1D>(const int i, const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mFwd.size());
      assert(this->mFwd.at(k).size() > 0);
      assert(i < this->mFwd.at(k).size());

      return this->mFwd.at(k)(i);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATF1D>(const int i) const
   {
      return this->idx<Dimensions::Data::DATF1D>(i, 0);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATB1D>(const int i, const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mBwd.size());
      assert(this->mBwd.at(k).size() > 0);
      assert(i < this->mBwd.at(k).size());

      return this->mBwd.at(k)(i);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATB1D>(const int i) const
   {
      return this->idx<Dimensions::Data::DATB1D>(i, 0);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DAT2D>(const int i, const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mIdx2D.size());
      assert(this->mIdx2D.at(k).size() > 0);
      assert(i < this->mIdx2D.at(k).size());

      return this->mIdx2D.at(k)(i);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DAT2D>(const int i) const
   {
      return this->idx<Dimensions::Data::DAT2D>(i, 0);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DAT3D>(const int i) const
   {
      // Check for correct sizes
      assert(this->mIdx3D.size() > 0);
      assert(i < this->mIdx3D.size());

      return this->mIdx3D(i);
   }

   ArrayI TransformResolution::mode(const int i) const
   {
      ArrayI mode(4); // 0 -> index 3D, 1 -> index 2D, 2 -> mode 3D, 3 -> mode 2D

      int current = 0;
      for(mode(0) = 0; mode(0) < this->mIdx3D.size(); ++mode(0))
      {
         mode(1) = i-current;
         if(mode(1) < this->mIdx2D.at(mode(0)).size())
         {
            mode(2) = this->mIdx3D(mode(0));
            mode(3) = this->mIdx2D.at(mode(0))(mode(1));
            current = i;
            break;
         } else
         {
            current += this->mIdx2D.at(mode(0)).size();
         }
      }

      // Safety assert
      assert(current == i);

      return mode;
   }

   void TransformResolution::clearIndexes()
   {
      // Create forward indexes
      std::vector<ArrayI>().swap(this->mFwd);

      // Create forward indexes
      std::vector<ArrayI>().swap(this->mBwd);

      // Clear indexes of second dimension
      std::vector<ArrayI>().swap(this->mIdx2D);

      // Clear indexes of third dimension
      this->mIdx3D.resize(0);
   }

}
