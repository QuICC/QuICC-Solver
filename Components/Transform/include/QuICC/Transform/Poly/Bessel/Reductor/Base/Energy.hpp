/**
 * @file Energy.hpp
 * @brief Implementation of the Bessel based energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_ENERGY_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_ENERGY_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/EnergyReductor.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/Power.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   template <class Impl>
   class Energy;

   /**
    * @brief Implementation of the Bessel based energy operator
    */
   template <>
   class Energy<base_t>: public EnergyReductor<Power<base_t>>
   {
      public:
         /**
          * @brief Constructor
          */
         Energy();

         /**
          * @brief Destructor
          */
         virtual ~Energy() = default;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_ENERGY_HPP
