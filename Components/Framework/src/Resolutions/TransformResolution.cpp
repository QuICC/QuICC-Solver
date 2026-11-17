/**
 * @file TransformResolution.cpp
 * @brief Source of the resolution object for a transform
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Resolutions/TransformResolution.hpp"

namespace QuICC {

   TransformResolution::TransformResolution(const std::vector<int>& globalDims, const std::vector<std::vector<std::vector<int> > >& fwd, const std::vector<std::vector<std::vector<int> > >& bwd, const std::vector<std::vector<int> >& idx2D, const std::vector<int>& idx3D)
      : mGlobalDims(globalDims), mFwd(fwd), mBwd(bwd), mIdx2D(idx2D), mIdx3D(idx3D), mDimF1D(), mDimB1D(), mDim2D(), mDim3D(idx3D.size())
   {
      // Initialise the dimensions
      this->initDimensions();
   }

   bool TransformResolution::isCleared() const
   {
      return (this->mIdx3D.size() == 0);
   }

   void TransformResolution::initDimensions()
   {
      for(int k = 0; k < this->mDim3D; k++)
      {
         this->mDim2D.push_back(this->mIdx2D.at(k).size());
      }

      this->mDimF1D.reserve(this->mDim3D);
      this->mDimB1D.reserve(this->mDim3D);
      for(int k = 0; k < this->mDim3D; k++)
      {
         this->mDimF1D.push_back(std::vector<int>());
         this->mDimB1D.push_back(std::vector<int>());
         for(int j = 0; j < this->mDim2D.at(k); j++)
         {
            this->mDimF1D.at(k).push_back(this->mFwd.at(k).at(j).size());
            this->mDimB1D.at(k).push_back(this->mBwd.at(k).at(j).size());
         }
      }
   }

   template <> int TransformResolution::dim<Dimensions::Data::DATF1D>(const int j, const int k) const
   {
      assert(this->mDimF1D.size() > static_cast<std::size_t>(k));
      assert(this->mDimF1D.at(k).size() > static_cast<std::size_t>(j));

      return this->mDimF1D.at(k).at(j);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DATF1D>() const
   {
      return this->dim<Dimensions::Data::DATF1D>(0,0);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DATB1D>(const int j, const int k) const
   {
      assert(this->mDimB1D.size() > static_cast<std::size_t>(k));
      assert(this->mDimB1D.at(k).size() > static_cast<std::size_t>(j));

      return this->mDimB1D.at(k).at(j);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DATB1D>() const
   {
      return this->dim<Dimensions::Data::DATB1D>(0,0);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DAT2D>(const int k) const
   {
      assert(this->mDim2D.size() > static_cast<std::size_t>(k));

      return this->mDim2D.at(k);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DAT2D>() const
   {
      return this->dim<Dimensions::Data::DAT2D>(0);
   }

   template <> int TransformResolution::dim<Dimensions::Data::DAT3D>() const
   {
      return this->mDim3D;
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATF1D>(const int i, const int j, const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mFwd.size());
      assert(this->mFwd.at(k).size() > 0);
      assert(static_cast<std::size_t>(i) < this->mFwd.at(k).at(j).size());
      assert(static_cast<std::size_t>(j) < this->mFwd.at(k).size());

      return this->mFwd.at(k).at(j).at(i);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATF1D>(const int i) const
   {
      return this->idx<Dimensions::Data::DATF1D>(i, 0, 0);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATB1D>(const int i, const int j, const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mBwd.size());
      assert(this->mBwd.at(k).size() > 0);
      assert(static_cast<std::size_t>(j) < this->mBwd.at(k).size());
      assert(static_cast<std::size_t>(i) < this->mBwd.at(k).at(j).size());

      return this->mBwd.at(k).at(j).at(i);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DATB1D>(const int i) const
   {
      return this->idx<Dimensions::Data::DATB1D>(i, 0, 0);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DAT2D>(const int i, const int k) const
   {
      // Check for correct sizes
      assert(k >= 0);
      assert(static_cast<size_t>(k) < this->mIdx2D.size());
      assert(this->mIdx2D.at(k).size() > 0);
      assert(static_cast<std::size_t>(i) < this->mIdx2D.at(k).size());

      return this->mIdx2D.at(k).at(i);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DAT2D>(const int i) const
   {
      return this->idx<Dimensions::Data::DAT2D>(i, 0);
   }

   template <> int TransformResolution::idx<Dimensions::Data::DAT3D>(const int i) const
   {
      // Check for correct sizes
      assert(this->mIdx3D.size() > 0);
      assert(static_cast<std::size_t>(i) < this->mIdx3D.size());

      return this->mIdx3D.at(i);
   }

   ArrayI TransformResolution::mode(const int i) const
   {
      ArrayI mode(4); // 0 -> index 3D, 1 -> index 2D, 2 -> mode 3D, 3 -> mode 2D

      int current = 0;
      for(mode(0) = 0; static_cast<std::size_t>(mode(0)) < this->mIdx3D.size(); ++mode(0))
      {
         mode(1) = i-current;
         if(static_cast<std::size_t>(mode(1)) < this->mIdx2D.at(mode(0)).size())
         {
            mode(2) = this->mIdx3D.at(mode(0));
            mode(3) = this->mIdx2D.at(mode(0)).at(mode(1));
            current = i;
            break;
         }
         else
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
      // Clear forward indexes
      std::vector<std::vector<std::vector<int> > >().swap(this->mFwd);

      // Clear forward indexes
      std::vector<std::vector<std::vector<int> > >().swap(this->mBwd);

      // Clear indexes of second dimension
      std::vector<std::vector<int> >().swap(this->mIdx2D);

      // Clear indexes of third dimension
      this->mIdx3D.resize(0);
   }

   std::shared_ptr<CscMetadata> TransformResolution::viewMeta(const Dimensions::Data::Id id) const
   {
      auto spMeta = std::make_shared<CscMetadata>();
      auto&& meta = *spMeta;

      meta.global2D = this->mGlobalDims.at(2);
      meta.global3D = this->mGlobalDims.at(3);

      auto it3D = this->mIdx3D.cbegin();
      auto itD2D = this->mDim2D.cbegin();
      std::pair<std::vector<std::vector<int>>::const_iterator, std::vector<std::vector<int>>::const_iterator> itD1D;
      if(id == Dimensions::Data::DATF1D)
      {
         meta.global1D = this->mGlobalDims.at(0);
         itD1D.first = this->mDimF1D.cbegin();
         itD1D.second = this->mDimF1D.cend();
      }
      else if(id == Dimensions::Data::DATB1D)
      {
         meta.global1D = this->mGlobalDims.at(1);
         itD1D.first = this->mDimB1D.cbegin();
         itD1D.second = this->mDimB1D.cend();
      }
      else
      {
         throw std::logic_error("Can only generate metadata for DATF1D and DATB1D");
      }

      meta.ptr2D.reserve(meta.global3D + 1);
      meta.ptr2D.push_back(0);
      std::uint32_t sze3D = 0;
      std::uint32_t sze2D = 0;
      for(std::uint32_t i = 0; i < meta.global3D; i++)
      {
         if(it3D != this->mIdx3D.end() && i == static_cast<uint32_t>(*it3D))
         {
            auto s = *itD2D;
            meta.ptr2D.push_back(meta.ptr2D.back() + s);

            it3D++;
            itD2D++;
            sze2D += s;
            sze3D++;
         }
         else
         {
            meta.ptr2D.push_back(meta.ptr2D.back());
         }
      }
      assert(sze3D == this->mIdx3D.size());

      // Store 1D dimensions
      meta.dim1D.reserve(sze2D);
      for(;itD1D.first != itD1D.second; itD1D.first++)
      {
         for(auto&& d: *(itD1D.first))
         {
            meta.dim1D.push_back(d);
         }
      }

      // Store 2D indices
      meta.idx2D.reserve(sze2D);
      for(auto&& d: this->mIdx2D)
      {
         for(auto&& i: d)
         {
            meta.idx2D.push_back(i);
         }
      }

      return spMeta;
   }

}
