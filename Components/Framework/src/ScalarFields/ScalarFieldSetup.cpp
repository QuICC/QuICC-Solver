/**
 * @file ScalarFieldSetup.cpp
 * @brief Source of the scalar field setup
 */

// Debug includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"

// Project includes
//
//
namespace QuICC {

namespace Datatypes {

   ScalarFieldSetup::ScalarFieldSetup(SharedArrayI spDim1D, SharedArrayI spDim2D, const int dim3D)
      : mspDim1D(spDim1D), mspDim2D(spDim2D), mDim3D(dim3D), mDataRows(0), mDataCols(0)
   {
      // Safety assertions
      assert(dim3D > 0 || (dim3D == 0 && this->mspDim2D->size() == 0 && this->mspDim1D->size() == 0));
      assert(this->mspDim2D->size() == this->mDim3D);
      assert(this->mspDim1D->size() == this->mDim3D);

      if(dim3D > 0)
      {
         assert(this->mspDim2D->minCoeff() > 0);
         assert(this->mspDim1D->minCoeff() > 0);

         this->mDataRows = this->mspDim1D->maxCoeff();
         this->mDataCols = this->mspDim2D->sum();
      }
   }

   ScalarFieldSetup::~ScalarFieldSetup()
   {
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

}
}
