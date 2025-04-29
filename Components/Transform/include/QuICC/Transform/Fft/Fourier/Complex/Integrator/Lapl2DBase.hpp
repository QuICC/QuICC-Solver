/**
 * @file Lapl2D.hpp
 * @brief Implementation of the Fourier based 2D Laplacian integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_LAPL2DBASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_LAPL2DBASE_HPP

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
   class Lapl2D;

   /**
    * @brief Implementation of the Fourier based 2D Laplacian integrator
    */
   template<>
   class Lapl2D<base_t>: public IComplexIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         Lapl2D() = default;

         /**
          * @brief Destructor
          */
         ~Lapl2D() = default;

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

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_LAPL2DBASE_HPP
