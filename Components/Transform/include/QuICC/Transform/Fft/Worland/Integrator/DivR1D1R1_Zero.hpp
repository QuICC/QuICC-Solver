/**
 * @file DivR1D1R1_Zero.hpp
 * @brief Implementation of the Worland based DivR1D1R1_Zero integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_DIVR1D1R1_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_DIVR1D1R1_ZERO_HPP

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
#include "QuICC/Transform/Fft/Worland/Integrator/DivR1D1R1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   /**
    * @brief Implementation of the Worland based DivR1D1R1_Zero integrator
    */
   class DivR1D1R1_Zero: public DivR1D1R1
   {
      public:
         /**
          * @brief Constructor
          */
         DivR1D1R1_Zero();

         /**
          * @brief Destructor
          */
         virtual ~DivR1D1R1_Zero();

      protected:
         /**
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const override;

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_DIVR1D1R1_ZERO_HPP
