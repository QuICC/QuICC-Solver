/**
 * @file PowerD1R1.hpp
 * @brief Implementation of the Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_POWERD1R1_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_POWERD1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tags.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/IEnergyWrapper.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/PowerD1R1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   template <class Impl>
   class PowerD1R1;

   /**
    * @brief Implementation of the Worland based R^2 energy operator
    */
   template <>
   class PowerD1R1<base_t>: public IEnergyWrapper<Poly::Worland::Reductor::PowerD1R1<Poly::Worland::base_t>>
   {
      public:
         /**
          * @brief Constructor
          */
         PowerD1R1() = default;

         /**
          * @brief Destructor
          */
         ~PowerD1R1() = default;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_POWERD1R1_HPP
