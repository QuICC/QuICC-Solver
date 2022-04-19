/** 
 * @file IMixedBackend.hpp
 * @brief Interface for a generic mixed cuFFT based backend 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_IMIXEDBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_IMIXEDBACKEND_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/ICuFftBackend.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Interface for a generic mixed cuFFT based backend
    */ 
   class IMixedBackend: public ICuFftBackend 
   {
      public:
         /// Typedef for the configuration class
         typedef Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         IMixedBackend();

         /**
          * @brief Constructor
          */
         IMixedBackend(const int nStreams, const int nBatches);

         /**
          * @brief Destructor
          */
         virtual ~IMixedBackend();
         
         /**
          * @brief Initialise the cuFFT transforms
          */
         virtual void init(const SetupType& setup) const;
         
      protected:
         /**
          * @brief Minimum number of blocks per batch
          */
         static const int MIN_BATCH_BLOCKSIZE;

         /**
          * @brief Get array of positive k
          */
         Array positiveK() const;

         /**
          * @brief Temporary GPU forward data
          */
         mutable std::vector<cufftDoubleReal*>  mcuFwd;

         /**
          * @brief Temporary GPU backward data
          */
         mutable std::vector<cufftDoubleComplex*>  mcuBwd;

         /**
          * @brief Padding size
          */
         mutable int mSpecSize;

         /**
          * @brief Forward FFT size
          */
         mutable int mFwdSize;

         /**
          * @brief Backward FFT size
          */
         mutable int mBwdSize;

         /**
          * @brief FFT blocksiz
          */
         mutable int mBlockSize;

         /**
          * @brief Temporary data
          */
         mutable MatrixZ  mTmp;

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_IMIXEDBACKEND_HPP
