/**
 * @file R1_Zero.hpp
 * @brief Implementation of the Worland based R1_Zero integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_R1_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_R1_ZERO_HPP

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
#include "QuICC/Transform/Fft/Worland/Integrator/R1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   /**
    * @brief Implementation of the Worland based R1_Zero integrator
    */
   class R1_Zero: public R1
   {
      public:
         /**
          * @brief Constructor
          */
         R1_Zero();

         /**
          * @brief Destructor
          */
         virtual ~R1_Zero();

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

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_R1_ZERO_HPP
