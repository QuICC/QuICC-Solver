/** 
 * @file PowerR2.hpp
 * @brief Implementation of the Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_POWERR2_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_POWERR2_HPP

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
#include "QuICC/Transform/Poly/Worland/Reductor/PowerR2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   /**
    * @brief Implementation of the Worland based R^2 energy operator
    */ 
   class PowerR2: public IEnergyWrapper<Poly::Worland::Reductor::PowerR2>
   {
      public:
         /**
          * @brief Constructor
          */
         PowerR2() = default;

         /**
          * @brief Destructor
          */
         virtual ~PowerR2() = default;
         
      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_POWERR2_HPP
