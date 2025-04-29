/**
 * @file P_Clean.hpp
 * @brief Implementation of the Fourier based P integrator, but 0 mode is cleaned
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_P_CLEANBASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_P_CLEANBASE_HPP

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
   class P_Clean;

   /**
    * @brief Implementation of the Fourier based P integrator, but 0 mode is cleaned
    */
   template <>
   class P_Clean<base_t>: public IComplexIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         P_Clean() = default;

         /**
          * @brief Destructor
          */
         ~P_Clean() = default;

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

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_P_CLEANBASE_HPP
