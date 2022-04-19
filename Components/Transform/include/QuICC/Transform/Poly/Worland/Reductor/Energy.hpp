/** 
 * @file Energy.hpp
 * @brief Implementation of the Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGY_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGY_HPP

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
#include "QuICC/Transform/Poly/Worland/Reductor/Power.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   /**
    * @brief Implementation of the Worland based energy operator
    */ 
   class Energy: public EnergyReductor<Power>
   {
      public:
         /**
          * @brief Constructor
          */
         Energy();

         /**
          * @brief Destructor
          */
         virtual ~Energy();
         
      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGY_HPP
