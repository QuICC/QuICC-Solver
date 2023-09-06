/**
 * @file P.hpp
 * @brief Implementation of the Fourier based P integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_PBASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_PBASE_HPP

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

   template<class Impl>
   class P;

   /**
    * @brief Implementation of the Fourier based P integrator
    */
   template<>
   class P<base_t>: public IMixedIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         P() = default;

         /**
          * @brief Destructor
          */
         ~P() = default;

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

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_PBASE_HPP
