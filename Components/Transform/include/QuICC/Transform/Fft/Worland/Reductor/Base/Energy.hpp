/**
 * @file Energy.hpp
 * @brief Implementation of the Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_ENERGY_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_ENERGY_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tags.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/IEnergyWrapper.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Base/Energy.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   template <class Impl>
   class Energy;

   /**
    * @brief Implementation of the Worland based R^2 energy operator
    */
   template <>
   class Energy<base_t>: public IEnergyWrapper<Poly::Worland::Reductor::Energy<Poly::Worland::base_t>>
   {
      public:
         /**
          * @brief Constructor
          */
         Energy() = default;

         /**
          * @brief Destructor
          */
         ~Energy() = default;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_BASE_ENERGY_HPP
