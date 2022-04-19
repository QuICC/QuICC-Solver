/**
 * @file Library.hpp
 * @brief Static interface to the global features of the CUDA cuFFT library  
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_LIBRARY_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_LIBRARY_HPP

// Configuration includes
//

// System includes
//

// External includes
//
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>

#ifdef QUICC_DEBUG
   #include <helper_cuda.h>
#else
   #define checkCudaErrors(x) x
#endif //QUICC_DEBUG

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Static interface to the global features of the CUDA cuFFT library 
    */ 
   class Library
   {
      public:
         /**
          * @brief Initialize the library
          */
         static void init();

         /**
          * @brief Register object using the library
          */
         static void registerUser();

         /**
          * @brief Unregister object using the library
          */
         static void unregisterUser();

         /**
          * @brief Cleanup the library
          */
         static void cleanup();

      private:
         /**
          * @brief Counter for the number of active objects
          */
         static int sCounter; 

         /**
          * @brief Empty constructor
          */
         Library();

         /**
          * @brief Empty destructor
          */
         ~Library(); 
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_LIBRARY_HPP
