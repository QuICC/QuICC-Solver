/**
 * @file IChebyshevBackend.hpp
 * @brief Interface for a generic Chebyshev FFTW based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_ICHEBYSHEVBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_ICHEBYSHEVBACKEND_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Backend/StorageKind.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/IFftwBackend.hpp"
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

/// This namespace provides all the Fftw related code.
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
          * @brief Copy input and pad
          *
          * @param tmp temporary storage
          * @param in input spectral coefficients
          */
         void input(Matrix& tmp, const Matrix& in) const;

         /**
          * @brief Copy input real or imaginary part and shift
          *
          * @param tmp temporary storage
          * @param in input spectral coefficients
          * @param shift
          */
         void input(Matrix& tmp, const Matrix& in, const int shift) const;

         /**
          * @brief Copy input real or imaginary part and pad
          *
          * @param tmp temporary storage
          * @param in input spectral coefficients
          * @param useReal flag to extract real or im part
          */
         void input(Matrix& tmp, const MatrixZ& in, const bool useReal) const;

         /**
          * @brief Copy input real or imaginary part and shift
          *
          * @param tmp temporary storage
          * @param in input spectral coefficients
          * @param shift
          * @param useReal flag to extract real or im part
          */
         void input(Matrix& tmp, const MatrixZ& in, const int shift, const bool useReal) const;

         /**
          * @brief Apply padding
          */
         virtual void applyPadding(Matrix& rData, const int extraRows = 0) const;

         /**
          * @brief Get the temporary storage
          *
          * @param getOut return input or ouput storage
          */
         virtual Matrix& getStorage(const StorageKind = StorageKind::in) const;

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
          * @brief Padding size
          */
         mutable int mPadSize = 0;

         /**
          * @brief Temporary data
          */
         mutable Matrix  mTmp;

         /**
          * @brief Temporary data for component wise operations
          */
         mutable Matrix  mTmpComp;

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_ICHEBYSHEVBACKEND_HPP
