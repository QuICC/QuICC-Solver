/** 
 * @file IChebyshevIntegrator.hpp
 * @brief Interface for a generic Chebyshev FFT based integrator 
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ICHEBYSHEVINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ICHEBYSHEVINTEGRATOR_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"
#include "QuICC/Transform/Fft/Backend/ChebyshevIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   /**
    * @brief Interface for a generic Chebyshev FFT based integrator
    */ 
   class IChebyshevIntegrator: public IChebyshevOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IChebyshevIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~IChebyshevIntegrator();

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
         Backend::ChebyshevIntegrator mBackend;

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyPreOperator(Matrix& rOut, const Matrix& in) const = 0;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(Matrix& rOut) const = 0;

         /**
          * @brief Apply pre FFT operator for component wise operations
          *
          * @param in   Input values
          */
         virtual void applyPreOperator(const MatrixZ& in, const bool useReal) const = 0;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(MatrixZ& rOut, const bool useReal) const = 0;

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

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ICHEBYSHEVINTEGRATOR_HPP
