/** 
 * @file P.hpp
 * @brief Implementation of the Fourier based mixed P projector (i.e. from phyiscal to modal space)
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_P_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_P_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Pect includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/IMixedProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   /**
    * @brief Implementation of the Fourier based mixed P projector
    * (i.e. from phyiscal to modal space)
    */ 
   class P: public IMixedProjector
   {
      public:
         /**
          * @brief Constructor
          */
         P();

         /**
          * @brief Destructor
          */
         ~P();

      protected:

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param in   Input mods values
          * @param out  Copied and padded input
          */
         void applyPreOperator(MatrixZ& out, const MatrixZ& in) const final;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_P_HPP
