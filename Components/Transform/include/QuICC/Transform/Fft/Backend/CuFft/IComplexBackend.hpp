/** 
 * @file IComplexBackend.hpp
 * @brief Backend for a generic complex cuFFT based backend 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_ICOMPLEXBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_ICOMPLEXBACKEND_HPP

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
#include "QuICC/Transform/Fft/Fourier/Complex/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Backend for a generic complex cuFFT based backend
    */ 
   class IComplexBackend: public ICuFftBackend
   {
      public:
         /// Typedef for the configuration class
         typedef Fourier::Complex::Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef Fourier::Complex::SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         IComplexBackend();

         /**
          * @brief Constructor
          */
         IComplexBackend(const int nStreams, const int nBatches);

         /**
          * @brief Destructor
          */
         virtual ~IComplexBackend();
         
         /**
          * @brief Initialise the cuFFT transforms
          */
         virtual void init(const SetupType& setup) const;

         /**
          * @brief Setup the mean blocks
          */
         void initMeanBlocks(const MatrixI& idBlocks) const;

         /**
          * @brief Set input data pointers for FFT (uses internal pointer for output)
          */
         virtual void input(const MHDComplex* in) const;

         /**
          * @brief Set output data pointers for FFT (uses internal pointer for input)
          */
         virtual void output(MHDComplex* out) const;

         /**
          * @brief Set input and output data pointers for FFT
          */
         virtual void io(MHDComplex* out, const MHDComplex* in) const;

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
          * @brief Get array of negative k
          */
         Array negativeK() const;

         /**
          * @brief Temporary GPU forward data
          */
         mutable std::vector<cufftDoubleComplex*>  mcuFwd;

         /**
          * @brief Temporary GPU backward data
          */
         mutable std::vector<cufftDoubleComplex*>  mcuBwd;

         /**
          * @brief Length of positive frequencies
          */
         mutable int mPosN;

         /**
          * @brief Length of negative frequencies
          */
         mutable int mNegN;

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
          * @brief Input data pointer
          */
         mutable const MHDComplex* mpIn;

         /**
          * @brief Out data pointer
          */
         mutable MHDComplex* mpOut;

         /**
          * @brief Temporary data
          */
         mutable MatrixZ  mTmp;

         /**
          * @brief Storage for the mean block sizes
          */
         mutable std::vector<std::pair<int,int> > mMeanBlocks;

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_ICOMPLEXBACKEND_HPP
