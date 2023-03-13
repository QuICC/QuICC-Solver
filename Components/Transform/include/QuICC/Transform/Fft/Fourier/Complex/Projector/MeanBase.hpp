/**
 * @file Mean.hpp
 * @brief Implementation of the Fourier based mean projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_MEANBASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_MEANBASE_HPP

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
   template<class>
   class Mean;

   /**
    * @brief Implementation of the Fourier based mean projector
    */
   template<>
   class Mean<base_t>: public IComplexProjector
   {
      public:
         /**
          * @brief Constructor
          */
         Mean() = default;

         /**
          * @brief Destructor
          */
         ~Mean() = default;

      protected:

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void applyPreOperator(MatrixZ& rOut, const MatrixZ& in) const;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_MEANBASE_HPP
