/**
 * @file P_Zero.hpp
 * @brief Implementation of the Worland based P_Zero integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_P_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_P_ZERO_HPP

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
   class P_Zero;

   /**
    * @brief Implementation of the Worland based P integrator and zero l = 0 mode
    */
   template <>
   class P_Zero<base_t>: public P<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         P_Zero();

         /**
          * @brief Destructor
          */
         ~P_Zero() = default;

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

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_P_ZERO_HPP
