/**
 * @file EnergyD1R1.hpp
 * @brief Implementation of the Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_ENERGYD1R1_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_ENERGYD1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tags.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/IEnergyWrapper.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Base/EnergyD1R1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   template <class Impl>
   class EnergyD1R1;

   /**
    * @brief Implementation of the Worland based R^2 energy operator
    */
   template <>
   class EnergyD1R1<base_t>: public IEnergyWrapper<Poly::Worland::Reductor::EnergyD1R1<Poly::Worland::base_t>>
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyD1R1() = default;

         /**
          * @brief Destructor
          */
         ~EnergyD1R1() = default;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_ENERGYD1R1_HPP
