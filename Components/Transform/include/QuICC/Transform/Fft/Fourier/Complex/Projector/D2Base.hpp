/**
 * @file D2.hpp
 * @brief Implementation of the Fourier based D^2 projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_D2BASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_D2BASE_HPP


// System includes
//


// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/IComplexProjector.hpp"
#include "QuICC/Transform/Fft/Fourier/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   template<class Impl>
   class D2;

   /**
    * @brief Implementation of the Fourier based D^2 projector
    */
   template<>
   class D2<base_t>: public IComplexProjector
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
          * @param rOut Output values
          * @param in   Input values
          */
         void applyPreOperator(MatrixZ& rOut, const MatrixZ& in) const final;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_D2BASE_HPP
