/** 
 * @file IChebyshevBackend.hpp
 * @brief Interface for a generic Chebyshev cuFFT based integrator 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_ICHEBYSHEVBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_ICHEBYSHEVBACKEND_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Interface for a generic Chebyshev cuFFT based integrator
    */ 
   class IChebyshevBackend: public ICuFftBackend
   {
      public:
         /// Typedef for the configuration class
         typedef Chebyshev::Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef Chebyshev::SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         IChebyshevBackend();

         /**
          * @brief Constructor
          */
         IChebyshevBackend(const int nStreams, const int nBatches);

         /**
          * @brief Destructor
          */
         virtual ~IChebyshevBackend();
         
         /**
          * @brief Initialise the FFT transforms
          */
         virtual void init(const SetupType& setup) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         virtual void io(double* out, const double* in) const;
         
      protected:
         /**
          * @brief Minimum number of blocks per batch
          */
         static const int MIN_BATCH_BLOCKSIZE;

         /**
          * @brief Temporary GPU forward data
          */
         mutable std::vector<cufftDoubleReal*>  mcuFwd;

         /**
          * @brief Temporary GPU backward data
          */
         mutable std::vector<cufftDoubleComplex*>  mcuBwd;

         /**
          * @brief Temporary GPU forward data
          */
         mutable std::vector<cufftDoubleReal*>  mcuWork;

         /**
          * @brief Spec size
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
          * @brief Padding size
          */
         mutable int mBlockSize;

         /**
          * @brief Temporary data
          */
         mutable Matrix  mTmp;

         /**
          * @brief Temporary data for component wise operations
          */
         mutable Matrix  mTmpComp;

         /**
          * @brief Input data pointer
          */
         mutable const double* mpIn;

         /**
          * @brief Out data pointer
          */
         mutable double* mpOut;

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_ICHEBYSHEVBACKEND_HPP
