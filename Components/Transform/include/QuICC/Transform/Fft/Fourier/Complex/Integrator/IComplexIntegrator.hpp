/**
 * @file IComplexIntegrator.hpp
 * @brief Interface for a generic complex FFTW based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_ICOMPLEXINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_ICOMPLEXINTEGRATOR_HPP

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
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/IComplexOperator.hpp"
#include "QuICC/Transform/Fft/Backend/ComplexIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   /**
    * @brief Interface for a generic complex FFTW based integrator
    */
   class IComplexIntegrator: public IComplexOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IComplexIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~IComplexIntegrator();

         /**
          * @brief Compute transform C2C
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const MatrixZ& in) const override;

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const override;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const override;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const override;

         /**
          * @brief dealias modes
          *
          * @param truncated  values
          * @param extended   values
          */
         void dealias(MatrixZ& deAliased, const MatrixZ& aliased) const final;


      protected:
         /**
          * @brief FFT backend
          */
         Backend::ComplexIntegrator mBackend;

      private:
         /**
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const override;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(MatrixZ& rOut) const = 0;

         /**
          * @brief Compute transform C2R (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute transform R2C (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const Matrix& in) const override;

         /**
          * @brief Compute transform R2R (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const Matrix& in) const override;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_ICOMPLEXINTEGRATOR_HPP
