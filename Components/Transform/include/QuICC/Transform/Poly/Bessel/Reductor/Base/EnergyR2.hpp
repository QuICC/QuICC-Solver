/**
 * @file EnergyR2.hpp
 * @brief Implementation of the Bessel based R^2 energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_ENERGYR2_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_ENERGYR2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/EnergyReductor.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/PowerR2.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   template <class Impl>
   class EnergyR2;

   /**
    * @brief Implementation of the Bessel based R^2 energy operator
    */
   template <>
   class EnergyR2<base_t>: public EnergyReductor<PowerR2<base_t>>
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyR2();

         /**
          * @brief Destructor
          */
         virtual ~EnergyR2() = default;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_ENERGYR2_HPP
