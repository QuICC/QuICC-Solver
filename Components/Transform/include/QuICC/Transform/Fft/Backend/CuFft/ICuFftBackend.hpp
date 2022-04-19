/** 
 * @file ICuFftBackend.hpp
 * @brief Interface for a generic cuFFT based backend 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_ICUFFTBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_ICUFFTBACKEND_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Backend/CuFft/Library.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Interface for a generic cuFFT backend
    */ 
   class ICuFftBackend
   {
      public:
         /**
          * @brief Constructor
          */
         ICuFftBackend();

         /**
          * @brief Constructor
          */
         ICuFftBackend(const int nStreams, const int nBatches);

         /**
          * @brief Destructor
          */
         virtual ~ICuFftBackend();

         /**
          * @brief Apply FFT
          */
         virtual void applyFft() const = 0;
         
      protected:
         /**
          * @brief Synchronize all streams
          */
         void synchronize() const;

         /**
          * @brief Synchronize single streams
          */
         void synchronize(const int id) const;

         /**
          * @brief Synchronize single streams
          */
         cudaStream_t stream(const int id) const;

         /**
          * @brief CUDA streams
          */
         mutable std::vector<cudaStream_t>  mStreams;

         /**
          * @brief Plans for the transform
          */
         mutable std::vector<cufftHandle>   mPlans;

         /**
          * @brief Number of CUDA streams
          */
         const int mNStreams;

         /**
          * @brief Number of subbatches
          */
         mutable int mNBatches;

      private:
         /**
          * @brief Initialise the cuFFT library
          */
         void initLibrary() const;

         /**
          * @brief Cleanup the FFT
          */
         void cleanupFft();

         /**
          * @brief Cleanup the cuFFT library on destruction
          */
         void cleanupLibrary();
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_ICUFFTBACKEND_HPP
