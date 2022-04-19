/**
 * @file Tools.hpp
 * @brief Definition of some useful constants and tools for FFT
 */

#ifndef QUICC_TRANSFORM_FFT_TOOLS_HPP
#define QUICC_TRANSFORM_FFT_TOOLS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/Tools.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

   /**
    * @brief Contains some useful constants and tools for FFTW
    */
   class Tools
   {
      public:
         /**
          * @brief Compute the dealiased size for FFT (R<->R, or Z<->Z)
          *
          * @param size Size to dealias
          */
         static int dealiasFft(const int size);

         /**
          * @brief Compute the dealiased size for FFT (Z<->Z)
          *
          * @param size Size to dealias
          */
         static int dealiasComplexFft(const int size);

         /**
          * @brief Compute the dealiased size for FFT (R<->Z, or Z<->R)
          *
          * @param size Size to dealias
          */
         static int dealiasMixedFft(const int size);

         /**
          * @brief Compute the dealiased size for cosine FFT
          *
          * @param size Size to dealias
          */
         static int dealiasCosFft(const int size);

         /**
          * @brief Optimise the FFT sizes
          *
          * @param size Current size of the FFT
          */
         static int optimizeFft(const int size);
         
      protected:
         /// Tools backend
         typedef Backend::Fftw::Tools Backend;

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

#endif // QUICC_TRANSFORM_FFT_TOOLS_HPP
