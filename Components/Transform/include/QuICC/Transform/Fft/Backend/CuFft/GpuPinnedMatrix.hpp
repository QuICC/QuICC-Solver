/**
 * @file GpuPinnedMatrix.hpp
 * @brief Simple struct to hold Matrix information in GPU pinned memory
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_GPUPINNNEDMATRIX_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_GPUPINNNEDMATRIX_HPP

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
    * @brief Simple struct to hold matrix information in GPU pinned memory
    */
   class GpuPinnedMatrix
   {
      public:
         /**
          * @brief Flag to keep size during reshape
          */
         static const int KEEP_SIZE = -1;

         /**
          * @brief Empty constructor
          */
         GpuPinnedMatrix();

         /**
          * @brief Allocate rows*cols GPU data
          */
         GpuPinnedMatrix(const int rows, const int cols);

         /**
          * @brief Copy constructor
          */
         GpuPinnedMatrix(const GpuPinnedMatrix& o);

         /**
          * @brief Move constructor
          */
         GpuPinnedMatrix(GpuPinnedMatrix&& o);

         /**
          * @brief Empty Destructor
          */
         ~GpuPinnedMatrix();

         /**
          * @brief Copy assigment
          */
         GpuPinnedMatrix& operator=(const GpuPinnedMatrix& o);

         /**
          * @brief Move assigment
          */
         GpuPinnedMatrix& operator=(GpuPinnedMatrix&& o);

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

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_GPUPINNNEDMATRIX_HPP
