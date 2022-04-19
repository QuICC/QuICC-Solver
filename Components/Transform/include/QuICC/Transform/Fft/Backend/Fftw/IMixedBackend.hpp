/**
 * @file IMixedBackend.hpp
 * @brief Interface for a generic mixed FFTW based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_IMIXEDBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_IMIXEDBACKEND_HPP

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
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Interface for a generic mixed FFTW based integrator
    */
   class IMixedBackend: public IFftwBackend
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
          * @brief Destructor
          */
         virtual ~IMixedBackend();

         /**
          * @brief Initialise the FFTW transforms
          */
         virtual void init(const SetupType& setup) const;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const;

      protected:
         /**
          * @brief Get array of positive k
          */
         Array positiveK() const;

         /**
          * @brief Padding size
          */
         mutable int mSpecSize;

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

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_IMIXEDBACKEND_HPP
