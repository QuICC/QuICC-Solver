/**
 * @file GpuCsrMatrix.hpp
 * @brief Simple struct to hold CSR Matrix information for GPU
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_GPUCSRMATRIX_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_GPUCSRMATRIX_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Simple struct to hold CSR matrix information
    */
   class GpuCsrMatrix
   {
      public:
         /**
          * @brief Flag to keep size during reshape
          */
         static const int KEEP_SIZE = -1;

         /**
          * @brief Empty constructor
          */
         GpuCsrMatrix();

         /**
          * @brief Allocate rows*cols GPU data
          */
         GpuCsrMatrix(const int rows, const int cols, const int nnz);

         /**
          * @brief Copy constructor
          */
         GpuCsrMatrix(const GpuCsrMatrix& o);

         /**
          * @brief Move constructor
          */
         GpuCsrMatrix(GpuCsrMatrix&& o);

         /**
          * @brief Empty Destructor
          */
         ~GpuCsrMatrix();

         /**
          * @brief Copy assigment
          */
         GpuCsrMatrix& operator=(const GpuCsrMatrix& o);

         /**
          * @brief Move assigment
          */
         GpuCsrMatrix& operator=(GpuCsrMatrix&& o);

         /**
          * @brief CSR data pointer
          */
         double* data();

         /**
          * @brief CSR data pointer
          */
         const double* data() const;

         /**
          * @brief CSR row pointer
          */
         int* rowPtr();

         /**
          * @brief CSR row pointer pointer
          */
         const int* rowPtr() const;

         /**
          * @brief CSR column index pointer
          */
         int* colIdx();

         /**
          * @brief CSR column index pointer
          */
         const int* colIdx() const;

         /**
          * @brief Number of rows
          */
         int rows() const;

         /**
          * @brief Number of rows
          */
         int cols() const;

         /**
          * @brief Number of nonzero entries
          */
         int nnz() const;

         /**
          * @brief Reshape if allocated memory is big enough
          */
         void reshape(int rows, int cols, int nnz);

         /**
          * @brief Allocate GPU memory
          */
         void allocate(const int rows, const int cols, const int nnz);

         /**
          * @brief Free GPU memory
          */
         void free();
         
      protected:

      private:
         /**
          * @brief Data pointer
          */
         double* mData;

         /**
          * @brief Row pointer pointer
          */
         int* mRowPtr;

         /**
          * @brief Column indexes pointer
          */
         int* mColIdx;

         /**
          * @brief Number of rows of allocated memory
          */
         int mGpuRows;

         /**
          * @brief Number of nonzero of allocated memory
          */
         int mGpuNnz;

         /**
          * @brief Number of rows
          */
         int mRows;

         /**
          * @brief Number of rows
          */
         int mCols;

         /**
          * @brief Number of nonzero entries
          */
         int mNnz;


   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_GPUCSRMATRIX_HPP
