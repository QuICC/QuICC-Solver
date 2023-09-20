/**
 * @file P.hpp
 * @brief Implementation of the Fourier based mixed P projector (i.e. from phyiscal to modal space)
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_PBASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_PBASE_HPP

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
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/IMixedProjector.hpp"
#include "QuICC/Transform/Fft/Fourier/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   template<class Impl>
   class P;

   /**
    * @brief Implementation of the Fourier based mixed P projector
    * (i.e. from modal to physical space)
    */
   template<>
   class P<base_t>: public IMixedProjector
   {
      public:
         /**
          * @brief Constructor
          */
         P() = default;

         /**
          * @brief Destructor
          */
         ~P() = default;

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
