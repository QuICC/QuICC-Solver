/** 
 * @file D1_Neg.hpp
 * @brief Implementation of the Fourier based D* integrator, but 0 mode is -P integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_D1_NEG_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_D1_NEG_HPP

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
    * @brief Implementation of the Fourier based D* integrator, but 0 mode is -P integrator
    */ 
   class D1_Neg: public IComplexIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         D1_Neg();

         /**
          * @brief Destructor
          */
         virtual ~D1_Neg();
         
      protected:

      private:
         /**
          * @brief Initialize operator with mean blocks
          */
         virtual void initOperator() const;

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

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_D1_NEG_HPP
