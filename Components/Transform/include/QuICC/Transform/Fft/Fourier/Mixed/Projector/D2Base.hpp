/** 
 * @file D2.hpp
 * @brief Implementation of the Fourier based D^2 projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D2BASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D2BASE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/IMixedProjector.hpp"
#include "QuICC/Transform/Fft/Fourier/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   template<class Impl>
   class D2;

   /**
    * @brief Implementation of the Fourier based D^2 projector
    */ 
   template<>
   class D2<base_t>: public IMixedProjector
   {
      public:
         /**
          * @brief Constructor
          */
         D2() = default;

         /**
          * @brief Destructor
          */
         ~D2() = default;
         
      protected:

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param in   Input values
          * @param out  Scaled input
          */
         void applyPreOperator(MatrixZ& out, const MatrixZ& in) const final;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D2BASE_HPP
