/** 
 * @file ChebyshevIntegrator.hpp
 * @brief Interface for a generic Chebyshev FFTW based integrator 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_CHEBYSHEVINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_CHEBYSHEVINTEGRATOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/IChebyshevBackend.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Interface for a generic Chebyshev FFTW based integrator
    */ 
   class ChebyshevIntegrator: public IChebyshevBackend
   {
      public:
         /**
          * @brief Constructor
          */
         ChebyshevIntegrator();

         /**
          * @brief Destructor
          */
         ~ChebyshevIntegrator();
         
         /**
          * @brief Initialise the FFTW transforms
          */
         void init(const SetupType& setup) const final;

         /**
          * @brief Set output
          */
         void output(Matrix& rOut) const;

         /**
          * @brief Set output
          */
         void output(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const;

         /**
          * @brief Set output multiplied aby spectral operator
          */
         void outputSpectral(Matrix& rOut) const;

         /**
          * @brief Set output multiplied aby spectral operator
          */
         void outputSpectral(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const;

         /**
          * @brief Set spectral operator
          */
         void setSpectralOperator(const SparseMatrix& mat) const;

         /**
          * @brief Set spectral operator
          */
         void setMeanOperator(const SparseMatrix& mat) const;
         
      protected:

      private:
         /**
          * @brief Bwd size
          */
         mutable int mBwdSize;

         /**
          * @brief FFT scaling factor
          */
         mutable MHDFloat mFftScaling;

         /**
          * @brief Spectral operator
          */
         mutable SparseMatrix mSpecOp;

         /**
          * @brief Spectral operator
          */
         mutable SparseMatrix mMeanOp;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_CHEBYSHEVINTEGRATOR_HPP
