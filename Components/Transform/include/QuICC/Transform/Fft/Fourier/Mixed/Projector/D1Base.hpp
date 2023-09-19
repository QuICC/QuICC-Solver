/**
 * @file D1.hpp
 * @brief Implementation of the Fourier based D projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D1BASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D1BASE_HPP

// System includes
//

// Project includes
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
   class D1;

   /**
    * @brief Implementation of the Fourier based D projector
    */
   template<>
   class D1<base_t>: public IMixedProjector
   {
      public:
         /**
          * @brief Constructor
          */
         D1() = default;

         /**
          * @brief Destructor
          */
         ~D1() = default;

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

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D1BASE_HPP
