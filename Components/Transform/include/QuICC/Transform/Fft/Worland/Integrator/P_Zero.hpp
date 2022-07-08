/**
 * @file P_Zero.hpp
 * @brief Implementation of the Worland based P_Zero integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_P_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_P_ZERO_HPP

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
#include "QuICC/Transform/Fft/Worland/Integrator/P.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   /**
    * @brief Implementation of the Worland based P integrator and zero l = 0 mode
    */
   class P_Zero: public P
   {
      public:
         /**
          * @brief Constructor
          */
         P_Zero();

         /**
          * @brief Destructor
          */
         virtual ~P_Zero();

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

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_P_ZERO_HPP
