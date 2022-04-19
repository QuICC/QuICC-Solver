/**
 * @file GpuCsrMatrix.cu
 * @brief Simple CSR Matrix representation for GPU data
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
#include "QuICC/Transform/Fft/Backend/CuFft/GpuCsrMatrix.hpp"

// Project includes
//
#include "QuICC/Transform/Fft/Backend/CuFft/CheckCuda.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   GpuCsrMatrix::GpuCsrMatrix()
      : mData(nullptr), mRowPtr(nullptr), mColIdx(nullptr), mGpuRows(0), mGpuNnz(0), mRows(0), mCols(0), mNnz(0)
   {
   }

   GpuCsrMatrix::GpuCsrMatrix(const int rows, const int cols, const int nnz)
      : mData(nullptr), mRowPtr(nullptr), mColIdx(nullptr), mGpuRows(0), mGpuNnz(0), mRows(rows), mCols(cols), mNnz(nnz)
   {
      this->allocate(this->mRows, this->mCols, this->mNnz);
   }

   GpuCsrMatrix::GpuCsrMatrix(const GpuCsrMatrix& o)
      : mData(nullptr), mRowPtr(nullptr), mColIdx(nullptr), mGpuRows(o.mGpuRows), mGpuNnz(o.mGpuNnz), mRows(o.mRows), mCols(o.mCols), mNnz(o.mNnz)
   {
      this->allocate(this->mRows, this->mCols, this->mNnz);
   }

   GpuCsrMatrix::GpuCsrMatrix(GpuCsrMatrix&& o)
      : mData(std::exchange(o.mData,nullptr)), mRowPtr(std::exchange(o.mRowPtr,nullptr)), mColIdx(std::exchange(o.mColIdx,nullptr)), mGpuRows(std::exchange(o.mGpuRows,0)), mGpuNnz(std::exchange(o.mGpuNnz,0)), mRows(std::exchange(o.mRows,0)), mCols(std::exchange(o.mCols,0)), mNnz(std::exchange(o.mNnz,0))
   {
      this->allocate(this->mRows, this->mCols, this->mNnz);
   }

   GpuCsrMatrix& GpuCsrMatrix::operator=(const GpuCsrMatrix& o)
   {
      if(this != &o)
      {
         if(this->rows() != o.rows() || this->nnz() != o.nnz())
         {
            this->allocate(o.mRows, o.mCols, o.mNnz);
         } else
         {
            if(this->mCols != o.mCols)
            {
               this->mCols = o.mCols;
            }
         }
         CheckCuda(cudaMemcpy(this->mData, o.mData, this->mNnz*sizeof(double), cudaMemcpyDeviceToDevice), __LINE__);
         CheckCuda(cudaMemcpy(this->mRowPtr, o.mColIdx, (this->mRows+1)*sizeof(int), cudaMemcpyDeviceToDevice), __LINE__);
         CheckCuda(cudaMemcpy(this->mColIdx, o.mColIdx, this->mNnz*sizeof(int), cudaMemcpyDeviceToDevice), __LINE__);
      }

      return *this;
   }

   GpuCsrMatrix& GpuCsrMatrix::operator=(GpuCsrMatrix&& o)
   {
      if(this != &o)
      {
         this->free();
         this->mData = std::exchange(o.mData, nullptr);
         this->mRowPtr = std::exchange(o.mRowPtr, nullptr);
         this->mColIdx = std::exchange(o.mColIdx, nullptr);
         this->mGpuRows = std::exchange(o.mGpuRows, 0);
         this->mGpuNnz = std::exchange(o.mGpuNnz, 0);
         this->mRows = std::exchange(o.mRows, 0);
         this->mCols = std::exchange(o.mCols, 0);
         this->mNnz = std::exchange(o.mNnz, 0);
      }

      return *this;
   }

   GpuCsrMatrix::~GpuCsrMatrix()
   {
      this->free();
   }

   double* GpuCsrMatrix::data()
   {
      return this->mData;
   }

   const double* GpuCsrMatrix::data() const
   {
      return this->mData;
   }

   int* GpuCsrMatrix::rowPtr()
   {
      return this->mRowPtr;
   }

   const int* GpuCsrMatrix::rowPtr() const
   {
      return this->mRowPtr;
   }

   int* GpuCsrMatrix::colIdx()
   {
      return this->mColIdx;
   }

   const int* GpuCsrMatrix::colIdx() const
   {
      return this->mColIdx;
   }

   int GpuCsrMatrix::rows() const
   {
      return this->mRows;
   }

   int GpuCsrMatrix::cols() const
   {
      return this->mCols;
   }

   int GpuCsrMatrix::nnz() const
   {
      return this->mNnz;
   }

   void GpuCsrMatrix::reshape(int rows, int cols, int nnz)
   {
      if(rows == GpuCsrMatrix::KEEP_SIZE)
      {
         rows = this->mRows;
      }

      if(cols == GpuCsrMatrix::KEEP_SIZE)
      {
         cols = this->mCols;
      }

      if(nnz == GpuCsrMatrix::KEEP_SIZE)
      {
         nnz = this->mNnz;
      }

      if(this->mGpuRows >= rows && this->mGpuNnz >= nnz)
      {
         this->mRows = rows;
         this->mCols = cols;
         this->mNnz = nnz;
      } else
      {
         throw std::logic_error("Tried to resize GpuCsrMatrix with too little memory allocated!");
      }
   }

   void GpuCsrMatrix::allocate(const int rows, const int cols, const int nnz)
   {
      this->free();
      this->mGpuRows = rows;
      this->mGpuNnz = nnz;
      CheckCuda(cudaMalloc((void**)&(this->mData), sizeof(double)*this->mGpuNnz), __LINE__);
      CheckCuda(cudaMalloc((void**)&(this->mRowPtr), sizeof(int)*(this->mGpuRows+1)), __LINE__);
      CheckCuda(cudaMalloc((void**)&(this->mColIdx), sizeof(int)*this->mGpuNnz), __LINE__);
      this->mRows = rows;
      this->mCols = cols;
      this->mNnz = nnz;
   }

   void GpuCsrMatrix::free()
   {
      if(this->mData)
      {
         CheckCuda(cudaFree(this->mData), __LINE__);
      }
      if(this->mRowPtr)
      {
         CheckCuda(cudaFree(this->mRowPtr), __LINE__);
      }
      if(this->mColIdx)
      {
         CheckCuda(cudaFree(this->mColIdx), __LINE__);
      }
   }

}
}
}
}
}
