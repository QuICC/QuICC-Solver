/**
 * @file Power.hpp
 * @brief Implementation of the Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_POWER_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_POWER_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tags.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/IEnergyWrapper.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Base/Power.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   template <class Impl>
   class Power;

   /**
    * @brief Implementation of the Worland based R^2 energy operator
    */
   template <>
   class Power<base_t>: public IEnergyWrapper<Poly::Worland::Reductor::Power<Poly::Worland::base_t>>
   {
      public:
         /**
          * @brief Constructor
          */
         Power() = default;

         /**
          * @brief Destructor
          */
         ~Power() = default;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_POWER_HPP
