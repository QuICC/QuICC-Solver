/** 
 * @file D2.hpp
 * @brief Implementation of the Chebyshev based D^2 projector, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_D2_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_D2_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/IChebyshevProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   /**
    * @brief Implementation of the Chebyshev based D^2 projector, with linear map y = ax + b
    */ 
   class D2: public IChebyshevProjector
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
         /**
          * @brief Initialize storage
          */
         virtual void initBackend() const override;

      private:
         /**
          * @brief Initialize solver operators
          */
         virtual void initOperator() const override;

         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyPreOperator(Matrix& rOut, const Matrix& in) const override;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(Matrix& rOut) const override;

         /**
          * @brief Apply pre FFT operator for component wise openerations
          *
          * @param in   Input values
          * @param useReal Real vs Imag flag
          */
         virtual void applyPreOperator(const MatrixZ& in, const bool useReal) const override;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          * @param useReal Real vs Imag flag
          */
         virtual void applyPostOperator(MatrixZ& rOut, const bool useReal) const override;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_D2_HPP
