/**
 * @file Ds1Lapl2D.hpp
 * @brief Implementation of the Fourier based D(slow) of 2D laplacian projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_DS1LAPL2DBASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_DS1LAPL2DBASE_HPP

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
   class Ds1Lapl2D;

   /**
    * @brief Implementation of the Fourier based D(slow) of 2D laplacian projector
    */
   template<>
   class Ds1Lapl2D<base_t>: public IComplexProjector
   {
      public:
         /**
          * @brief Constructor
          */
         Ds1Lapl2D() = default;

         /**
          * @brief Destructor
          */
         ~Ds1Lapl2D() = default;

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

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_DS1LAPL2DBASE_HPP
