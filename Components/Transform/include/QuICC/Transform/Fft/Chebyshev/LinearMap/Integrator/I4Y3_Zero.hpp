/**
 * @file I4Y3_Zero.hpp
 * @brief Implementation of the Chebyshev based I^4 Y^3 integrator, but 0 mode is zeroed, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I4Y3_ZERO_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I4Y3_ZERO_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/IChebyshevIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   /**
    * @brief Implementation of the Chebyshev based I^4 Y^4 of 1/Y integrator, but 0 mode is zeroed, with linear map y = ax + b
    */
   class I4Y3_Zero: public IChebyshevIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         I4Y3_Zero();

         /**
          * @brief Destructor
          */
         ~I4Y3_Zero();

      protected:
         /**
          * @brief Sparse matrix operator
          */
         mutable SparseMatrix mOp;

      private:
         /**
          * @brief Initialize solver operators
          */
         void initOperator() const final;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         void applyPostOperator(Matrix& rOut) const final;

         /**
          * @brief Apply pre FFT operator for component wise openerations
          *
          * @param in   Input values
          * @param useReal Real vs Imag flag
          */
         void applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const final;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          * @param useReal Real vs Imag flag
          */
         void applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const final;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I4Y3_ZERO_HPP
