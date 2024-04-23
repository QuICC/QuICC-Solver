/**
 * @file EnergySLaplR2.hpp
 * @brief Implementation of the Bessel based spherical Laplacian R^2 energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_ENERGYSLAPLR2_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_ENERGYSLAPLR2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/EnergyReductor.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/PowerSLaplR2.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   template <class Impl>
   class EnergySLaplR2;

   /**
    * @brief Implementation of the Bessel based Spherical Laplacian R^2 energy operator
    */
   template <>
   class EnergySLaplR2<base_t>: public EnergyReductor<PowerSLaplR2<base_t>>
   {
      public:
         /**
          * @brief Constructor
          */
         EnergySLaplR2();

         /**
          * @brief Destructor
          */
         virtual ~EnergySLaplR2() = default;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_ENERGYSLAPLR2_HPP
