/** 
 * @file EnergyR2.hpp
 * @brief Implementation of the Worland based R^2 energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYR2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYR2_HPP

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
#include "QuICC/Transform/Poly/Worland/Reductor/EnergyReductor.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/PowerR2.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   /**
    * @brief Implementation of the Worland based R^2 energy operator
    */ 
   class EnergyR2: public EnergyReductor<PowerR2>
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyR2();

         /**
          * @brief Destructor
          */
         virtual ~EnergyR2();
         
      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYR2_HPP
