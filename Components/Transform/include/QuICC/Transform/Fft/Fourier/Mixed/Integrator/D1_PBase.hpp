/**
 * @file D1_PBase.hpp
 * @brief Implementation of the Fourier based D integrator, but 0 mode is P integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_PBASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_PBASE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/IMixedIntegrator.hpp"
#include "QuICC/Transform/Fft/Fourier/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Integrator {

   template<class>
   class D1_P;

   /**
    * @brief Implementation of the Fourier based D integrator but 0 mode is P integrator
    */
   template<>
   class D1_P<base_t>: public IMixedIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         D1_P() = default;

         /**
          * @brief Destructor
          */
         ~D1_P() = default;

      protected:

      private:
         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         void applyPostOperator(MatrixZ& rOut) const final;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_PBASE_HPP
