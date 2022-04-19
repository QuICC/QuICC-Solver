/** 
 * @file D3.hpp
 * @brief Implementation of the Fourier based D^3 projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D3_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D3_HPP

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
    * @brief Implementation of the Fourier based D^3 projector
    */ 
   class D3: public IMixedProjector
   {
      public:
         /**
          * @brief Constructor
          */
         D3();

         /**
          * @brief Destructor
          */
         virtual ~D3();
         
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

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D3_HPP
