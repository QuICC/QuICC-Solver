/** 
 * @file EnergyD1.hpp
 * @brief Implementation of the Chebyshev based D^1 energy reductor, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGYD1_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGYD1_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/IChebyshevEnergy.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   /**
    * @brief Implementation of the Chebyshev based D^1 energy reductor, with linear map y = ax + b
    */ 
   class EnergyD1: public IChebyshevEnergy
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyD1();

         /**
          * @brief Destructor
          */
         ~EnergyD1();
         
      protected:
         /**
          * @brief Initialise operator
          */
         void initOperator() const final;

         /**
          * @brief Initialize storage
          */
         void initBackend() const final;

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param in   Input values
          */
         void applyPreOperator(Matrix& tmp, const Matrix& in) const final;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         void applyPostOperator(Matrix& rOut, const Matrix& tmp) const final;

         /**
          * @brief Apply pre FFT operator for component wise openerations
          *
          * @param in   Input values
          * @param useReal Real vs Imag flag
          */
         void applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const final;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGYD1_HPP
