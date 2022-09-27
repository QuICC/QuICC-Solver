/** 
 * @file IComplexBackend.hpp
 * @brief Backend for a generic complex FFTW based projector 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_ICOMPLEXBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_ICOMPLEXBACKEND_HPP

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
#include "QuICC/Transform/Fft/Backend/Fftw/IFftwBackend.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Backend for a generic complex FFTW based projector
    */ 
   class IComplexBackend: public IFftwBackend
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
          * @brief Destructor
          */
         virtual ~IComplexBackend();
         
         /**
          * @brief Initialise the FFTW transforms
          */
         virtual void init(const SetupType& setup) const;

         /**
          * @brief Setup the mean blocks
          */
         void initMeanBlocks(const MatrixI& idBlocks) const;

         /**
          * @brief Get the temporary storage
          */
         MatrixZ& getStorage() const;

      protected:
         /**
          * @brief Get array of positive k
          */
         Array positiveK() const;

         /**
          * @brief Get array of negative k
          */
         Array negativeK() const;

         /**
          * @brief Length of positive frequencies
          */
         mutable int mPosN;

         /**
          * @brief Length of negative frequencies
          */
         mutable int mNegN;

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

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_ICOMPLEXBACKEND_HPP
