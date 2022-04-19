/**
 * @file Tools.hpp
 * @brief Definition of some useful constants and tools for CUDA cuFFT 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_TOOLS_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_TOOLS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Contains some useful constants and tools for CUDA cuFFT
    */
   class Tools
   {
      public:
         /**
          * @brief Optimise the FFT sizes
          *
          * @param size Current size of the FFT
          */
         static bool optimizeFft(const int size, int& opt);

         /**
          * @brief Standard dealiasing factor (usually 3/2)
          */
         static const double STD_DEALIASING;

         /**
          * @brief The real <-> complex fast Fourrier transform dealiasing factor
          */
         static const double MIXED_DEALIASING;

         /**
          * @brief Cosine dealiasing factor (usually 3/2)
          */
         static const double COS_DEALIASING;

         /**
          * @brief Maximul extension width to consider for optimization
          */
         static const double OPTIMIZATION_WIDTH;
         
      protected:

      private:

         /**
          * @brief Empty constructor
          */
         Tools();

         /**
          * @brief Empty Destructor
          */
         ~Tools();

   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_TOOLS_HPP
