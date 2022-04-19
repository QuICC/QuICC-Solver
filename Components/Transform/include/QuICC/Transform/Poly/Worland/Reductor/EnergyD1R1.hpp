/** 
 * @file EnergyD1R1.hpp
 * @brief Implementation of the Worland based D R energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYD1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYD1R1_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/EnergyReductor.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/PowerD1R1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   /**
    * @brief Implementation of the Worland based D R energy operator
    */ 
   class EnergyD1R1: public EnergyReductor<PowerD1R1>
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyD1R1();

         /**
          * @brief Destructor
          */
         virtual ~EnergyD1R1();
         
      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYD1R1_HPP
