/**
 * @file GpuPinnedMatrix.cu
 * @brief Simple Matrix representation for GPU data
 */

// System includes
//
#include <stdexcept>
#include <utility>

// External includes
//
#include <cublas_v2.h>

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/GpuPinnedMatrix.hpp"

// Project includes
//
#include "QuICC/Transform/Fft/Backend/CuFft/CheckCuda.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   GpuPinnedMatrix::GpuPinnedMatrix()
      : mData(nullptr), mMemSize(0), mRows(0), mCols(0)
   {
   }

   GpuPinnedMatrix::GpuPinnedMatrix(const int rows, const int cols)
      : mData(nullptr), mMemSize(0), mRows(rows), mCols(cols)
   {
      this->allocate(this->mRows, this->mCols);
   }

   GpuPinnedMatrix::GpuPinnedMatrix(const GpuPinnedMatrix& o)
      : mData(nullptr), mMemSize(o.mMemSize), mRows(o.mRows), mCols(o.mCols)
   {
      this->allocate(this->mRows, this->mCols);
   }

   GpuPinnedMatrix::GpuPinnedMatrix(GpuPinnedMatrix&& o)
      : mData(std::exchange(o.mData,nullptr)), mMemSize(std::exchange(o.mMemSize,0)), mRows(std::exchange(o.mRows,0)), mCols(std::exchange(o.mCols,0))
   {
      this->allocate(this->mRows, this->mCols);
   }

   GpuPinnedMatrix& GpuPinnedMatrix::operator=(const GpuPinnedMatrix& o)
   {
      if(this != &o)
      {
         if(this->size() != o.size())
         {
            this->allocate(o.mRows, o.mCols);
         } else
         {
            if(this->mRows != o.mRows)
            {
               this->mRows = o.mRows;
               this->mCols = o.mCols;
            }
         }
         cublasHandle_t handle;
         CheckCuda(cublasCreate(&handle), __LINE__);
         CheckCuda(cublasDcopy(handle, this->size(), o.mData, 1, this->mData, 1), __LINE__);
         CheckCuda(cublasDestroy(handle), __LINE__);
      }

      return *this;
   }

   GpuPinnedMatrix& GpuPinnedMatrix::operator=(GpuPinnedMatrix&& o)
   {
      if(this != &o)
      {
         this->free();
         this->mData = std::exchange(o.mData, nullptr);
         this->mMemSize = std::exchange(o.mMemSize, 0);
         this->mRows = std::exchange(o.mRows, 0);
         this->mCols = std::exchange(o.mCols, 0);
      }

      return *this;
   }

   GpuPinnedMatrix::~GpuPinnedMatrix()
   {
      this->free();
   }

   double* GpuPinnedMatrix::data()
   {
      return this->mData;
   }

   const double* GpuPinnedMatrix::data() const
   {
      return this->mData;
   }

   int GpuPinnedMatrix::rows() const
   {
      return this->mRows;
   }

   int GpuPinnedMatrix::cols() const
   {
      return this->mCols;
   }

   int GpuPinnedMatrix::size() const
   {
      return this->mRows*this->mCols;
   }

   void GpuPinnedMatrix::reshape(int rows, int cols)
   {
      if(rows == GpuPinnedMatrix::KEEP_SIZE)
      {
         rows = this->mRows;
      }

      if(cols == GpuPinnedMatrix::KEEP_SIZE)
      {
         cols = this->mCols;
      }

      if(this->mMemSize >= rows*cols)
      {
         this->mRows = rows;
         this->mCols = cols;
      } else
      {
         throw std::logic_error("Tried to resize GpuPinnedMatrix with too little memory allocated!");
      }
   }

   void GpuPinnedMatrix::resize(int rows, int cols)
   {
      if(rows == GpuPinnedMatrix::KEEP_SIZE)
      {
         rows = this->mRows;
      }

      if(cols == GpuPinnedMatrix::KEEP_SIZE)
      {
         cols = this->mCols;
      }

      this->mRows = rows;
      this->mCols = cols;

      if(this->mMemSize < rows*cols)
      {
         this->allocate(rows,cols);
      }
   }

   void GpuPinnedMatrix::allocate(const int rows, const int cols)
   {
      this->free();
      this->mMemSize = rows*cols;
      CheckCuda(cudaMallocHost((void**)&(this->mData), sizeof(double)*this->mMemSize), __LINE__);
      this->mRows = rows;
      this->mCols = cols;
   }

   void GpuPinnedMatrix::free()
   {
      if(this->mData)
      {
         CheckCuda(cudaFreeHost(this->mData), __LINE__);
      }
   }

}
}
}
}
}
