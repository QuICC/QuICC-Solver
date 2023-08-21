/**
 * @file DivR1_Zero.hpp
 * @brief Implementation of the Worland based 1/R projector and zero for l = 0
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_BASE_DIVR1_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_BASE_DIVR1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Tags.hpp"
#include "QuICC/Transform/Fft/Worland/Projector/IWorlandProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   template <class Impl>
   class DivR1_Zero;

   /**
    * @brief Implementation of the Worland based 1/R projector and zero for l = 0
    */
   template <>
   class DivR1_Zero<base_t>: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         DivR1_Zero();

         /**
          * @brief Destructor
          */
         ~DivR1_Zero() = default;

      protected:
         /**
          * @brief Initialise FFT backend
          */
         void initBackend() const final;

         /**
          * @brief Initialise FFT backend
          */
         void computeWorlandExpansion(const bool isEven) const final;

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

#endif // QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_BASE_DIVR1_ZERO_HPP
