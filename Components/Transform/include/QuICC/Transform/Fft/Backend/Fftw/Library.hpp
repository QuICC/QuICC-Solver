/**
 * @file Library.hpp
 * @brief Static interface to the global features of the FFTW library  
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_LIBRARY_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_LIBRARY_HPP

// Configuration includes
//

// System includes
//

// External includes
//
#include <fftw3.h>

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Static interface to the global features of the FFTW library 
    */ 
   class Library
   {
      public:
         /**
          * @brief Get the plan flag
          */
         static unsigned int planFlag();

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
          * @brief FFTW3 flags for the plan setup
          */
         static unsigned int  sPlanFlag;

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

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_LIBRARY_HPP
