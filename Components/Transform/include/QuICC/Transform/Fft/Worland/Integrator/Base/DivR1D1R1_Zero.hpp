/**
 * @file DivR1D1R1_Zero.hpp
 * @brief Implementation of the Worland based DivR1D1R1_Zero integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_DIVR1D1R1_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_DIVR1D1R1_ZERO_HPP

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
   class DivR1D1R1_Zero;

   /**
    * @brief Implementation of the Worland based DivR1D1R1_Zero integrator
    */
   template <>
   class DivR1D1R1_Zero<base_t>: public DivR1D1R1<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         DivR1D1R1_Zero();

         /**
          * @brief Destructor
          */
         ~DivR1D1R1_Zero() = default;

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

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_BASE_DIVR1D1R1_ZERO_HPP
