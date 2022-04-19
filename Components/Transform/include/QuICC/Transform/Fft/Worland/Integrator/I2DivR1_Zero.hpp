/**
 * @file I2DivR1_Zero.hpp
 * @brief Implementation of the Worland based I2DivR1 integrator and zero l = 0 mode
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_I2DIVR1_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_I2DIVR1_ZERO_HPP

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
#include "QuICC/Transform/Fft/Worland/Integrator/DivR1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   /**
    * @brief Implementation of the Worland based I2DivR1 integrator and zero l = 0 mode
    */
   class I2DivR1_Zero: public DivR1
   {
      public:
         /**
          * @brief Constructor
          */
         I2DivR1_Zero();

         /**
          * @brief Destructor
          */
         virtual ~I2DivR1_Zero();

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
         virtual void applyPostOperator(Matrix& rOut, const bool isEven) const override;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          * @param useReal Real vs Imag flag
          */
         virtual void applyPostOperator(MatrixZ& rOut, const bool useReal, const bool isEven) const override;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_I2DIVR1_ZERO_HPP
