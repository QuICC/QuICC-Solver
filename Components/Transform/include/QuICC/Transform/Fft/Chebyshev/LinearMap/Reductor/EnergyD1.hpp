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
         virtual ~EnergyD1();
         
      protected:
         /**
          * @brief Initialise operator
          */
         virtual void initOperator() const override;

         /**
          * @brief Initialize storage
          */
         virtual void initBackend() const override;

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param in   Input values
          */
         virtual void applyPreOperator(const Matrix& in) const override;

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
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGYD1_HPP
