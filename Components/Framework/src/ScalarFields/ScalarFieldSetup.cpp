/**
 * @file ScalarFieldSetup.cpp
 * @brief Source of the scalar field setup
 */

// System includes
//

// Project includes
//
#include "Memory/MemoryResource.hpp"
#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"

namespace QuICC {

namespace Datatypes {

   ScalarFieldSetup::ScalarFieldSetup(std::shared_ptr<CscMetadata> spMeta, std::shared_ptr<Memory::memory_resource> mem)
      : mspDim1D(0), mspDim2D(0), mDim3D(0), mDataRows(0), mDataCols(0), mspMeta(spMeta), mMem(mem)
   {
      // Count layers
      std::vector<int> layers;
      std::vector<int> dim1D;
      std::vector<int> dim2D;
      for(std::size_t i = 0; i < spMeta->ptr2D.size() - 1; i++)
      {
         if(spMeta->ptr2D.at(i+1) > spMeta->ptr2D.at(i))
         {
            layers.push_back(i);
            dim2D.push_back(spMeta->ptr2D.at(i+1) - spMeta->ptr2D.at(i));
            int max1D = 0;
            for(std::uint32_t j = spMeta->ptr2D.at(i); j < spMeta->ptr2D.at(i+1); j++)
            {
               max1D = std::max(max1D, static_cast<int>(spMeta->dim1D.at(j)));
            }
            dim1D.push_back(max1D);
         }
      }

      this->mDim3D = layers.size();
      assert(dim2D.size() == static_cast<std::size_t>(this->mDim3D));
      assert(dim1D.size() == static_cast<std::size_t>(this->mDim3D));

      mspDim1D = std::make_shared<ArrayI>(this->mDim3D);
      for(std::size_t i = 0; i < dim1D.size(); i++)
      {
         (*this->mspDim1D)(i) = dim1D.at(i);
      }
      mspDim2D = std::make_shared<ArrayI>(this->mDim3D);
      for(std::size_t i = 0; i < dim2D.size(); i++)
      {
         (*this->mspDim2D)(i) = dim2D.at(i);
      }

      // Safety assertions
      assert(this->mDim3D > 0 || (this->mDim3D == 0 && this->mspDim2D->size() == 0 && this->mspDim1D->size() == 0));
      assert(this->mspDim2D->size() == this->mDim3D);
      assert(this->mspDim1D->size() == this->mDim3D);

      if(this->mDim3D > 0)
      {
         assert(this->mspDim2D->minCoeff() > 0);
         assert(this->mspDim1D->minCoeff() > 0);

         this->mDataRows = this->mspDim1D->maxCoeff();
         this->mDataCols = this->mspDim2D->sum();
      }
   }

   int ScalarFieldSetup::dataRows() const
   {
      return this->mDataRows;
   }

   int ScalarFieldSetup::dataCols() const
   {
      return this->mDataCols;
   }

   int ScalarFieldSetup::colIdx(const int j, const int k) const
   {
      assert(this->mspDim2D->size() >= k);
      assert((*this->mspDim2D)(k) > j);

      return this->mspDim2D->head(k).sum() + j;
   }

   int ScalarFieldSetup::blockIdx(const int k) const
   {
      assert(this->mspDim2D->size() >= k);

      return this->mspDim2D->head(k).sum();
   }

   int ScalarFieldSetup::blockRows(const int k) const
   {
      assert(this->mspDim1D->size() > k);

      return (*this->mspDim1D)(k);
   }

   int ScalarFieldSetup::blockCols(const int k) const
   {
      assert(this->mspDim2D->size() > k);

      return (*this->mspDim2D)(k);
   }

   int ScalarFieldSetup::nBlock() const
   {
      return this->mDim3D;
   }

   std::shared_ptr<CscMetadata> ScalarFieldSetup::viewMeta() const
   {
      return this->mspMeta;
   }

   std::shared_ptr<Memory::memory_resource> ScalarFieldSetup::mem() const
   {
      return this->mMem;
   }

}
}
