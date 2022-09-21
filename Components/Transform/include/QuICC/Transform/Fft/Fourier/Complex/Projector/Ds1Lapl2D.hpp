/** 
 * @file Ds1Lapl2D.hpp
 * @brief Implementation of the Fourier based D(slow) of 2D laplacian projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_DS1LAPL2D_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_DS1LAPL2D_HPP

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
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/IComplexProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   /**
    * @brief Implementation of the Fourier based D(slow) of 2D laplacian projector
    */ 
   class Ds1Lapl2D: public IComplexProjector
   {
      public:
         /**
          * @brief Constructor
          */
         Ds1Lapl2D();

         /**
          * @brief Destructor
          */
         ~Ds1Lapl2D();
         
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

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_DS1LAPL2D_HPP
