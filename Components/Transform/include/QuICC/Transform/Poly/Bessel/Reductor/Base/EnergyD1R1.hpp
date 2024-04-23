/**
 * @file EnergyD1R1.hpp
 * @brief Implementation of the Bessel based D R energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_ENERGYD1R1_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_ENERGYD1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/EnergyReductor.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/PowerD1R1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   template <class Impl>
   class EnergyD1R1;


   /**
    * @brief Implementation of the Bessel based D R energy operator
    */
   template <>
   class EnergyD1R1<base_t>: public EnergyReductor<PowerD1R1<base_t>>
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyD1R1();

         /**
          * @brief Destructor
          */
         virtual ~EnergyD1R1() = default;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_ENERGYD1R1_HPP
