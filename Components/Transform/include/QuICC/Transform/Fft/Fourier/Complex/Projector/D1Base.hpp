/** 
 * @file D1.hpp
 * @brief Implementation of the Fourier based D projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_D1BASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_D1BASE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/IComplexProjector.hpp"
#include "QuICC/Transform/Fft/Fourier/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   template<class Impl>
   class D1;

   /**
    * @brief Implementation of the Fourier based D projector
    */ 
   template<>
   class D1<base_t>: public IComplexProjector
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

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_D1BASE_HPP
