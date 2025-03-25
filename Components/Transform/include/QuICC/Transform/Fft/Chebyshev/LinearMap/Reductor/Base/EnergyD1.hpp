/**
 * @file EnergyD1.hpp
 * @brief Implementation of the Chebyshev based D^1 energy reductor, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_BASE_ENERGYD1_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_BASE_ENERGYD1_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Tags.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/ILinearMapEnergy.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   template <typename Impl> class EnergyD1;

   /**
    * @brief Implementation of the Chebyshev based D^1 energy reductor, with linear map y = ax + b
    */
   template <>
   class EnergyD1<base_t>: public ILinearMapEnergy
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyD1() = default;

         /**
          * @brief Destructor
          */
         ~EnergyD1() = default;

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

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_BASE_ENERGYD1_HPP
