/**
 * @file P.hpp
 * @brief Implementation of the Chebyshev based P projector, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_P_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_P_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/IChebyshevProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   /**
    * @brief Implementation of the Chebyshev based P projector, with linear map y = ax + b
    */
   class P: public IChebyshevProjector
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
          * @param tmp Temporary padded modal values
          * @param in   Input values
          */
         void applyPreOperator(Matrix& tmp, const Matrix& in) const final;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         void applyPostOperator(Matrix& rOut) const final;

         /**
          * @brief Apply pre FFT operator for component wise openerations
          *
          * @param tmp Temporary padded modal values, either real or im part only
          * @param in Input values
          * @param useReal Real vs Imag flag
          */
         void applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const final;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Scaled values
          * @param tmp Real or imag coefficients
          * @param useReal Real vs Imag flag
          */
         void applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const final;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_P_HPP
