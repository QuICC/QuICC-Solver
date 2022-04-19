/** 
 * @file ChebyshevIntegrator.hpp
 * @brief Interface for a generic Chebyshev cuFFT based integrator 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_CHEBYSHEVINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_CHEBYSHEVINTEGRATOR_HPP

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
#include "QuICC/Transform/Fft/Backend/CuFft/IChebyshevBackend.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Interface for a generic Chebyshev cuFFT based integrator
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
         virtual ~ChebyshevIntegrator();
         
         /**
          * @brief Initialise the FFT transforms
          */
         virtual void init(const SetupType& setup) const override;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         virtual void io(double* out, const double* in) const override;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         virtual void io(Matrix& out, const Matrix& in) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         virtual void input(const MatrixZ& out, const bool useReal) const;

         /**
          * @brief Set output
          */
         void output(Matrix& rOut) const;

         /**
          * @brief Set output
          */
         void output(MatrixZ& rOut, const bool useReal) const;

         /**
          * @brief Set output multiplied aby spectral operator
          */
         void outputSpectral(Matrix& rOut) const;

         /**
          * @brief Set output multiplied aby spectral operator
          */
         void outputSpectral(MatrixZ& rOut, const bool useReal) const;

         /**
          * @brief Apply FFT
          */
         virtual void applyFft() const override;

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
          * @brief FFT scaling factor
          */
         mutable double mFftScaling;

         /**
          * @brief Map for output data
          */
         mutable Eigen::Map<Matrix> mOutMap;

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

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_CHEBYSHEVINTEGRATOR_HPP
