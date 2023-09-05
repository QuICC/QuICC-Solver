/**
 * @file PowerSLaplR2.hpp
 * @brief Implementation of the Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_POWERSLAPLR2_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_POWERSLAPLR2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tags.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/IEnergyWrapper.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/PowerSLaplR2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   template <class Impl>
   class PowerSLaplR2;

   /**
    * @brief Implementation of the Worland based R^2 energy operator
    */
   template <>
   class PowerSLaplR2<base_t>: public IEnergyWrapper<Poly::Worland::Reductor::PowerSLaplR2<Poly::Worland::base_t>>
   {
      public:
         /**
          * @brief Constructor
          */
         PowerSLaplR2() = default;

         /**
          * @brief Destructor
          */
         ~PowerSLaplR2() = default;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_POWERSLAPLR2_HPP
