/**
 * @file EnergySLaplR2.hpp
 * @brief Implementation of the Worland based spherical Laplacian R^2 energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_ENERGYSLAPLR2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_ENERGYSLAPLR2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Tags.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Base/EnergyReductor.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/PowerSLaplR2.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   template <class Impl>
   class EnergySLaplR2;

   /**
    * @brief Implementation of the Worland based Spherical Laplacian R^2 energy operator
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

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_ENERGYSLAPLR2_HPP
