/**
 * @file D1_Neg.hpp
 * @brief Implementation of the Fourier based D* integrator, but 0 mode is -P integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_D1_NEGBASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_D1_NEGBASE_HPP

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
   class D1_Neg;

   /**
    * @brief Implementation of the Fourier based D* integrator, but 0 mode is -P integrator
    */
   template<>
   class D1_Neg<base_t>: public IComplexIntegrator
   {
      public:

      protected:

      private:
         /**
          * @brief Initialize operator with mean blocks
          */
         void initOperator() const final;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         void applyPostOperator(MatrixZ& rOut) const final;

         /**
          * @brief Storage for the mean block sizes
          */
         mutable std::vector<std::pair<int,int> > mMeanBlocks;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_D1_NEGBASE_HPP
