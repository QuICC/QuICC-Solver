/** 
 * @file Energy.hpp
 * @brief Implementation of the Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_ENERGY_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_ENERGY_HPP

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
#include "QuICC/Transform/Fft/Worland/Reductor/IEnergyWrapper.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Energy.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   /**
    * @brief Implementation of the Worland based R^2 energy operator
    */ 
   class Energy: public IEnergyWrapper<Poly::Worland::Reductor::Energy>
   {
      public:
         /**
          * @brief Constructor
          */
         Energy() = default;

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

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_ENERGY_HPP
