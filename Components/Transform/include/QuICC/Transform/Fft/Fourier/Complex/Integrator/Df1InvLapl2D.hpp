/** 
 * @file Df1InvLapl2D.hpp
 * @brief Implementation of the Fourier based  D(fast) of inverse 2D Laplacian integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_DF1INVLAPL2D_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_DF1INVLAPL2D_HPP

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
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/IComplexIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   /**
    * @brief Implementation of the Fourier based D(fast) of inverse 2D Laplacian integrator
    */ 
   class Df1InvLapl2D: public IComplexIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         Df1InvLapl2D();

         /**
          * @brief Destructor
          */
         virtual ~Df1InvLapl2D();
         
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

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_DF1INVLAPL2D_HPP
