/**
 * @file Df1InvLapl2D.hpp
 * @brief Implementation of the Fourier based  D(fast) of inverse 2D Laplacian integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_DF1INVLAPL2DBASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_DF1INVLAPL2DBASE_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/IComplexIntegrator.hpp"
#include "QuICC/Transform/Fft/Fourier/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {
   template<class>
   class Df1InvLapl2D;

   /**
    * @brief Implementation of the Fourier based D(fast) of inverse 2D Laplacian integrator
    */
   template<>
   class Df1InvLapl2D<base_t>: public IComplexIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         Df1InvLapl2D() = default;

         /**
          * @brief Destructor
          */
         ~Df1InvLapl2D() = default;

      protected:

      private:
         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         void applyPostOperator(MatrixZ& rOut) const final;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_DF1INVLAPL2DBASE_HPP
