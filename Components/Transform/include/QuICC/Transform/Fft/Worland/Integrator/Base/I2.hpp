/**
 * @file I2.hpp
 * @brief Implementation of the Worland based I2 integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_I2_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_I2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tags.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/P.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   template <class Impl>
   class I2;

   /**
    * @brief Implementation of the Worland based I2 integrator
    */
   template <>
   class I2<base_t>: public P<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         I2();

         /**
          * @brief Destructor
          */
         virtual ~I2() = default;

      protected:
         /**
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const override;

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

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_I2_Zero_HPP
