/**
 * @file IMixedIntegrator.hpp
 * @brief Interface for a generic mixed FFTW based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_IMIXEDINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_IMIXEDINTEGRATOR_HPP

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
#include "QuICC/Transform/Fft/Fourier/Mixed/IMixedOperator.hpp"
#include "QuICC/Transform/Fft/Backend/MixedIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Integrator {

   /**
    * @brief Interface for a generic mixed FFTW based integrator
    */
   class IMixedIntegrator: public IMixedOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IMixedIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~IMixedIntegrator();

         /**
          * @brief Compute transform R2C
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const Matrix& in) const override;

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
         virtual void dealias(MatrixZ& deAliased, const MatrixZ& aliased) const override;

      protected:
         /**
          * @brief FFT backend
          */
         Backend::MixedIntegrator mBackend;

         /**
          * @brief Compute transform C2C or R2R componentwise (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute transform C2R (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute transform R2R (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const Matrix& in) const override;

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
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_IMIXEDINTEGRATOR_HPP
