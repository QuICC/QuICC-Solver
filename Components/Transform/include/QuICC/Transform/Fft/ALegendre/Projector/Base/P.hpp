/**
 * @file P.hpp
 * @brief Implementation of the ALegendre based P projector
 */

#ifndef QUICC_TRANSFORM_FFT_ALEGENDRE_PROJECTOR_BASE_P_HPP
#define QUICC_TRANSFORM_FFT_ALEGENDRE_PROJECTOR_BASE_P_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/ALegendre/Tags.hpp"
#include "QuICC/Transform/Fft/ALegendre/Projector/IALegendreProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace ALegendre {

namespace Projector {

   template <class Impl>
   class P;

   /**
    * @brief Implementation of the ALegendre based P projector
    */
   template <>
   class P<base_t>: public IALegendreProjector
   {
      public:
         /**
          * @brief Constructor
          */
         P();

         /**
          * @brief Destructor
          */
         ~P() = default;

      protected:
         /**
          * @brief Initialise FFT backend
          */
         void computeALegendreExpansion(const bool isEven) const final;

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void applyPreOperator(const Matrix& in, const bool isEven) const final;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         void applyPostOperator(Matrix& rOut, const bool isEven) const final;

         /**
          * @brief Apply pre FFT operator for component wise openerations
          *
          * @param in   Input values
          * @param useReal Real vs Imag flag
          */
         void applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const final;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          * @param useReal Real vs Imag flag
          */
         void applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const final;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_ALegendre_PROJECTOR_BASE_P_HPP
