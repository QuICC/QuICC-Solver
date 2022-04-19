/** 
 * @file I4Y3D1Y1_Zero.hpp
 * @brief Implementation of the Chebyshev based I^4 Y^4 of 1/Y D Y integrator, but 0 mode is zeroed, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I4Y3D1Y1_ZERO_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I4Y3D1Y1_ZERO_HPP

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
    * @brief Implementation of the Chebyshev based I^4 Y^4 of 1/Y D Y integrator, but 0 mode is zeroed, with linera map y = ax + b
    */ 
   class I4Y3D1Y1_Zero: public IChebyshevIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         I4Y3D1Y1_Zero();

         /**
          * @brief Destructor
          */
         virtual ~I4Y3D1Y1_Zero();
         
      protected:
         /**
          * @brief Sparse matrix operator
          */
         mutable SparseMatrix mOp;

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

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I4Y3D1Y1_ZERO_HPP
