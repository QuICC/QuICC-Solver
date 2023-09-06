/**
 * @file I2_Zero.hpp
 * @brief Implementation of the Worland based I2_Zero integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_I2_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_I2_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tags.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/I2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   template <class Impl>
   class I2_Zero;

   /**
    * @brief Implementation of the Worland based I2 integrator and zero l = 0 mode
    */
   template <>
   class I2_Zero<base_t>: public I2<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         I2_Zero();

         /**
          * @brief Destructor
          */
         ~I2_Zero() = default;

      protected:
         /**
          * @brief Initialise FFT backend
          */
         void initBackend() const final;

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_I2_ZERO_HPP
