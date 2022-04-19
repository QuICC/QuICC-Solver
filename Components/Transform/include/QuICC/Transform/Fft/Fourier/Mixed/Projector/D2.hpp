/** 
 * @file D2.hpp
 * @brief Implementation of the Fourier based D^2 projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D2_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D2_HPP

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
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/IMixedProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   /**
    * @brief Implementation of the Fourier based D^2 projector
    */ 
   class D2: public IMixedProjector
   {
      public:
         /**
          * @brief Constructor
          */
         D2();

         /**
          * @brief Destructor
          */
         virtual ~D2();
         
      protected:

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyPreOperator(Matrix& rOut, const MatrixZ& in) const;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(Matrix& rOut) const;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D2_HPP
