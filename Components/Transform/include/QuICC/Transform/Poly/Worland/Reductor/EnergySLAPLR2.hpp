/** 
 * @file EnergySLAPLR2.hpp
 * @brief Implementation of the Worland based spherical Laplacian R^2 energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYSLAPLR2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYSLAPLR2_HPP

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
#include "QuICC/Transform/Poly/Worland/Reductor/PowerSLAPLR2.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   /**
    * @brief Implementation of the Worland based Spherical Laplacian R^2 energy operator
    */ 
   class EnergySLAPLR2: public EnergyReductor<PowerSLAPLR2>
   {
      public:
         /**
          * @brief Constructor
          */
         EnergySLAPLR2();

         /**
          * @brief Destructor
          */
         virtual ~EnergySLAPLR2();
         
      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYSLAPLR2_HPP
