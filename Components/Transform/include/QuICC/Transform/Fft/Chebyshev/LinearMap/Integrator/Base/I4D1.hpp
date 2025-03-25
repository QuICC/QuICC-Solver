/**
 * @file I4D1.hpp
 * @brief Implementation of the Chebyshev based I^4 of D integrator, with y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_BASE_I4D1_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_BASE_I4D1_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Tags.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/ILinearMapIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   template <typename Impl> class I4D1;

   /**
    * @brief Implementation of the Chebyshev based I^4 of D integrator, with y = ax + b
    */
   template <>
   class I4D1<base_t> final : public ILinearMapIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         I4D1() = default;

         /**
          * @brief Destructor
          */
         ~I4D1() final = default;

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

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_BASE_I4D1_HPP
