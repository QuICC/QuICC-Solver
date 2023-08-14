/**
 * @file PowerD1R1.hpp
 * @brief Implementation of the Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_POWERD1R1_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_POWERD1R1_HPP

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
#include "QuICC/Transform/Poly/Worland/Reductor/PowerD1R1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   /**
    * @brief Implementation of the Worland based R^2 energy operator
    */
   class PowerD1R1: public IEnergyWrapper<Poly::Worland::Reductor::PowerD1R1<Poly::Worland::base_t>>
   {
      public:
         /**
          * @brief Constructor
          */
         PowerD1R1() = default;

         /**
          * @brief Destructor
          */
         virtual ~PowerD1R1() = default;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_POWERD1R1_HPP
