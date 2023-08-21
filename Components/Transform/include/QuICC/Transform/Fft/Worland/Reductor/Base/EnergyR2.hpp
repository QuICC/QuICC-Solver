/**
 * @file EnergyR2.hpp
 * @brief Implementation of the Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_ENERGYR2_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_ENERGYR2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tags.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/IEnergyWrapper.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Base/EnergyR2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   template <class Impl>
   class EnergyR2;

   /**
    * @brief Implementation of the Worland based R^2 energy operator
    */
   template <>
   class EnergyR2<base_t>: public IEnergyWrapper<Poly::Worland::Reductor::EnergyR2<Poly::Worland::base_t>>
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyR2() = default;

         /**
          * @brief Destructor
          */
         ~EnergyR2() = default;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_ENERGYR2_HPP
