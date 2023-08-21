/**
 * @file I4DivR1D1R1_Zero.hpp
 * @brief Implementation of the Worland based I4DivR1D1R1 integrator and zero l = 0 mode
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_I4DIVR1D1R1_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_I4DIVR1D1R1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tags.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/DivR1D1R1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   template <class Impl>
   class I4DivR1D1R1_Zero;

   /**
    * @brief Implementation of the Worland based I4DivR1D1R1 integrator and zero l = 0 mode
    */
   template <>
   class I4DivR1D1R1_Zero<base_t>: public DivR1D1R1<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         I4DivR1D1R1_Zero();

         /**
          * @brief Destructor
          */
         ~I4DivR1D1R1_Zero() = default;

      protected:
         /**
          * @brief Initialise FFT backend
          */
         void initBackend() const final;

      private:
         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         void applyPostOperator(Matrix& rOut, const bool isEven) const final;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          * @param useReal Real vs Imag flag
          */
         void applyPostOperator(MatrixZ& rOut, const bool useReal, const bool isEven) const final;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_I4DIVR1D1R1_ZERO_HPP
