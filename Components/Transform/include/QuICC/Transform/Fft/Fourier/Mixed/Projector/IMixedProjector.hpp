/**
 * @file IMixedProjector.hpp
 * @brief Interface for a generic mixed FFTW based projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_IMIXEDPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_IMIXEDPROJECTOR_HPP

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
#include "QuICC/Transform/Fft/Backend/MixedProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   /**
    * @brief Interface for a generic mixed FFTW based projector
    */
   class IMixedProjector: public IMixedOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IMixedProjector();

         /**
          * @brief Destructor
          */
         virtual ~IMixedProjector();

         /**
          * @brief Compute transform C2R
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const MatrixZ& in) const override;

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

      protected:
         /**
          * @brief FFT backend
          */
         Backend::MixedProjector mBackend;

      private:
         /**
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const override;

         /**
          * @brief Apply pre FFT operator
          *
          * @param out  Copied or scaled input
          * @param in   Input values
          */
         virtual void applyPreOperator(MatrixZ& out, const MatrixZ& in) const = 0;

         /**
          * @brief Compute transform C2C or R2R componentwise (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const MatrixZ& in) const override;

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

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_IMIXEDPROJECTOR_HPP
