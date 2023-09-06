/**
 * @file D2Base.hpp
 * @brief Implementation of the Fourier based D^2* integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D2BASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D2BASE_HPP

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
   class D2;

   /**
    * @brief Implementation of the Fourier based D^2* integrator
    */
   template<>
   class D2<base_t>: public IMixedIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         D2() = default;

         /**
          * @brief Destructor
          */
         ~D2() = default;

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

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D2BASE_HPP
