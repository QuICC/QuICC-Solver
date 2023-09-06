/** 
 * @file P.hpp
 * @brief Implementation of the Fourier based P projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_PBASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_PBASE_HPP

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
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/IComplexProjector.hpp"
#include "QuICC/Transform/Fft/Fourier/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   template<class Impl>
   class P;

   /**
    * @brief Implementation of the Fourier based P projector
    */ 
   template<>
   class P<base_t>: public IComplexProjector
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

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_P_HPP
