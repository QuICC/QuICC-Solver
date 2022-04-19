/**
 * @file I2D1_I2.hpp
 * @brief Implementation of the Chebyshev based I^2 of D integrator, but 0 mode is I^2 of P integrator, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I2D1_I2_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I2D1_I2_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/IChebyshevIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   /**
    * @brief Implementation of the Chebyshev based I^2 of D integrator, but 0 mode is I^2 of P integrator, with linear map y = ax + b
    */
   class I2D1_I2: public IChebyshevIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         I2D1_I2();

         /**
          * @brief Destructor
          */
         virtual ~I2D1_I2();

      protected:
         /**
          * @brief Sparse matrix operator
          */
         mutable SparseMatrix mOp;

         /**
          * @brief Sparse matrix operator for mean
          */
         mutable SparseMatrix mMeanOp;

      private:
         /**
          * @brief Initialize solver operators
          */
         virtual void initOperator() const;

         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyPreOperator(Matrix& rOut, const Matrix& in) const;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(Matrix& rOut) const;

         /**
          * @brief Apply pre FFT operator for component wise openerations
          *
          * @param in   Input values
          * @param useReal Real vs Imag flag
          */
         virtual void applyPreOperator(const MatrixZ& in, const bool useReal) const;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          * @param useReal Real vs Imag flag
          */
         virtual void applyPostOperator(MatrixZ& rOut, const bool useReal) const;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I2D1_I2_HPP
