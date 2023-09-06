/**
 * @file EnergyD1R1.hpp
 * @brief Implementation of the Worland based D R energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_ENERGYD1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_ENERGYD1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Tags.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Base/EnergyReductor.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/PowerD1R1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   template <class Impl>
   class EnergyD1R1;


   /**
    * @brief Implementation of the Worland based D R energy operator
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

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_BASE_ENERGYD1R1_HPP
