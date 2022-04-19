/** 
 * @file P.hpp
 * @brief Implementation of the Fourier based P projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_P_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_P_HPP

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

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   /**
    * @brief Implementation of the Fourier based P projector
    */ 
   class P: public IComplexProjector
   {
      public:
         /**
          * @brief Constructor
          */
         P();

         /**
          * @brief Destructor
          */
         virtual ~P();
         
      protected:

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyPreOperator(MatrixZ& rOut, const MatrixZ& in) const;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(MatrixZ& rOut) const;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_P_HPP
