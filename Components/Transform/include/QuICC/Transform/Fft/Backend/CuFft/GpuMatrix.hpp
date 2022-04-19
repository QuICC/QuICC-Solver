/**
 * @file GpuMatrix.hpp
 * @brief Simple struct to hold Matrix information for GPU
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_GPUMATRIX_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_GPUMATRIX_HPP

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
    * @brief Simple struct to hold matrix information
    */
   class GpuMatrix
   {
      public:
         /**
          * @brief Flag to keep size during reshape
          */
         static const int KEEP_SIZE = -1;

         /**
          * @brief Empty constructor
          */
         GpuMatrix();

         /**
          * @brief Allocate rows*cols GPU data
          */
         GpuMatrix(const int rows, const int cols);

         /**
          * @brief Copy constructor
          */
         GpuMatrix(const GpuMatrix& o);

         /**
          * @brief Move constructor
          */
         GpuMatrix(GpuMatrix&& o);

         /**
          * @brief Empty Destructor
          */
         ~GpuMatrix();

         /**
          * @brief Copy assigment
          */
         GpuMatrix& operator=(const GpuMatrix& o);

         /**
          * @brief Move assigment
          */
         GpuMatrix& operator=(GpuMatrix&& o);

         /**
          * @brief Data pointer
          */
         double* data();

         /**
          * @brief Data pointer
          */
         const double* data() const;

         /**
          * @brief Number of rows
          */
         int rows() const;

         /**
          * @brief Number of rows
          */
         int cols() const;

         /**
          * @brief Total size
          */
         int size() const;

         /**
          * @brief Reshape if allocated memory is big enough
          */
         void reshape(int rows, int cols);

         /**
          * @brief Reshape if allocated memory is big enough or resize memory otherwise
          */
         void resize(int rows, int cols);

         /**
          * @brief Allocate GPU memory
          */
         void allocate(const int rows, const int cols);

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
          * @brief Size of allocated memory
          */
         int mMemSize;

         /**
          * @brief Number of rows
          */
         int mRows;

         /**
          * @brief Number of rows
          */
         int mCols;


   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_GPUMATRIX_HPP
