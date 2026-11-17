/**
 * @file D4.hpp
 * @brief Implementation of the Fourier based D^4 projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D4_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D4_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/IMixedProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   /**
    * @brief Implementation of the Fourier based D^4 projector
    */
   class D4: public IMixedProjector
   {
      public:
         /**
          * @brief Constructor
          */
         D4() = default;

         /**
          * @brief Destructor
          */
         ~D4() = default;

      protected:

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param in   Input values
          * @param out  Scaled input
          */
         void applyPreOperator(MatrixZ& tmp, const Eigen::Ref<const MatrixZ>& in) const final;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D4_HPP
