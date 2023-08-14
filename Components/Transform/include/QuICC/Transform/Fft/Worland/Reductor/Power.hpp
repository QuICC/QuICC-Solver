/**
 * @file Power.hpp
 * @brief Implementation of the Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_POWER_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_POWER_HPP

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
#include "QuICC/Transform/Fft/Worland/Reductor/IEnergyWrapper.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Base/Power.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   /**
    * @brief Implementation of the Worland based R^2 energy operator
    */
   class Power: public IEnergyWrapper<Poly::Worland::Reductor::Power<Poly::Worland::base_t>>
   {
      public:
         /**
          * @brief Constructor
          */
         Power() = default;

         /**
          * @brief Destructor
          */
         virtual ~Power() = default;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_POWER_HPP
