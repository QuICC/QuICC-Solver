/** 
 * @file D1_P.hpp
 * @brief Implementation of the Fourier based D integrator, but 0 mode is P integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_P_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_P_HPP

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
    * @brief Implementation of the Fourier based D integrator but 0 mode is P integrator
    */ 
   class D1_P: public IMixedIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         D1_P();

         /**
          * @brief Destructor
          */
         virtual ~D1_P();
         
      protected:

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyPreOperator(MatrixZ& rOut, const Matrix& in) const;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(MatrixZ& rOut) const;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_P_HPP
