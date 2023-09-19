/**
 * @file Energy.hpp
 * @brief Implementation of the Chebyshev based energy reductor, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGY_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGY_HPP

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
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/IChebyshevEnergy.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   /**
    * @brief Implementation of the Chebyshev based energy reductor, with linear map y = ax + b
    */
   class Energy: public IChebyshevEnergy
   {
      public:
         /**
          * @brief Constructor
          */
         Energy();

         /**
          * @brief Destructor
          */
         ~Energy();

      protected:

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void applyPreOperator(Matrix& tmp, const Matrix& in) const  final;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         void applyPostOperator(Matrix& rOut, const Matrix& tmp) const  final;

         /**
          * @brief Apply pre FFT operator for component wise openerations
          *
          * @param in   Input values
          * @param useReal Real vs Imag flag
          */
         void applyPreOperator(Matrix& tmp, const MatrixZ& in, const  bool useReal) const final;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGY_HPP
