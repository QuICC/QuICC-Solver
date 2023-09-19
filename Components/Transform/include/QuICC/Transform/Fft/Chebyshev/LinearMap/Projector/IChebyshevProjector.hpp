/**
 * @file IChebyshevProjector.hpp
 * @brief Interface for a generic Chebyshev FFT based projector
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ICHEBYSHEVPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ICHEBYSHEVPROJECTOR_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"
#include "QuICC/Transform/Fft/Backend/ChebyshevProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   /**
    * @brief Interface for a generic Chebyshev FFT based projector
    */
   class IChebyshevProjector: public IChebyshevOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IChebyshevProjector();

         /**
          * @brief Destructor
          */
         virtual ~IChebyshevProjector();

         /**
          * @brief Compute transform R2R componentwise
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute transform R2R
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const Matrix& in) const override;

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
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const override;

         /**
          * @brief FFT backend
          */
         Backend::ChebyshevProjector mBackend;

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param tmp  Temporary padded modal values
          * @param in   Input values
          */
         virtual void applyPreOperator(Matrix& tmp, const Matrix& in) const = 0;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(Matrix& rOut) const = 0;

         /**
          * @brief Apply pre FFT operator for component wise operations
          *
          * @param tmp Temporary padded modal values, either real or im part only
          * @param in Input values
          * @param useReal 1 -> extract real part, 0 -> extract im part
          */
         virtual void applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const = 0;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const = 0;

         /**
          * @brief Compute transform R2C (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const Matrix& in) const override;

         /**
          * @brief Compute transform C2R (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const MatrixZ& in) const override;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ICHEBYSHEVPROJECTOR_HPP
