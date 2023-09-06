/**
 * @file R1_Zero.hpp
 * @brief Implementation of the Worland based R1_Zero integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_R1_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_R1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tags.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/R1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   template <class Impl>
   class R1_Zero;

   /**
    * @brief Implementation of the Worland based R1_Zero integrator
    */
   template <>
   class R1_Zero<base_t>: public R1<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         R1_Zero();

         /**
          * @brief Destructor
          */
         ~R1_Zero() = default;

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

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_R1_ZERO_HPP
