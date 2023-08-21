/**
 * @file DivR1_Zero.hpp
 * @brief Implementation of the Worland based DivR1_Zero integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_DIVR1_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_DIVR1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tags.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/DivR1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   template <class Impl>
   class DivR1_Zero;

   /**
    * @brief Implementation of the Worland based DivR1_Zero integrator
    */
   template <>
   class DivR1_Zero<base_t>: public DivR1<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         DivR1_Zero();

         /**
          * @brief Destructor
          */
         ~DivR1_Zero() = default;

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

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_DIVR1_ZERO_HPP
