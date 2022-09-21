/** 
 * @file D1_Neg.hpp
 * @brief Implementation of the Fourier based D* integrator, but 0 mode is -P integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_NEG_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_NEG_HPP

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
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/IMixedIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Integrator {

   /**
    * @brief Implementation of the Fourier based D* integrator, but 0 mode is -P integrator
    */ 
   class D1_Neg: public IMixedIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         D1_Neg();

         /**
          * @brief Destructor
          */
         ~D1_Neg();
         
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

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_NEG_HPP
