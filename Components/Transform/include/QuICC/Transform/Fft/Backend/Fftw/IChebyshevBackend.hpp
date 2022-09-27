/** 
 * @file IChebyshevBackend.hpp
 * @brief Interface for a generic Chebyshev FFTW based integrator 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_ICHEBYSHEVBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_ICHEBYSHEVBACKEND_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Interface for a generic Chebyshev FFTW based integrator
    */ 
   class IChebyshevBackend: public IFftwBackend
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
          * @brief Destructor
          */
         virtual ~IChebyshevBackend();
         
         /**
          * @brief Initialise the FFTW transforms
          */
         virtual void init(const SetupType& setup) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         virtual void io(MHDFloat* out, const MHDFloat* in) const;

         /**
          * @brief Get the temporary storage
          *
          * @param getOut return input or ouput storage
          */
         Matrix& getStorage(const bool getOut = false) const;
         
      protected:
         /**
          * @brief Spec size
          */
         mutable int mSpecSize;

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
         mutable const MHDFloat* mpIn;

         /**
          * @brief Out data pointer
          */
         mutable MHDFloat* mpOut;

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_ICHEBYSHEVBACKEND_HPP
