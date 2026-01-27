/**
 * @file D1.hpp
 * @brief Implementation of the ALegendre based D1 integrator
 */

#ifndef QUICC_TRANSFORM_FFT_ALEGENDRE_INTEGRATOR_BASE_D1_HPP
#define QUICC_TRANSFORM_FFT_ALEGENDRE_INTEGRATOR_BASE_D1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/ALegendre/Tags.hpp"
#include "QuICC/Transform/Fft/ALegendre/Integrator/IALegendreIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class D1;

   /**
    * @brief Implementation of the ALegendre based D1 integrator
    */
   template <>
   class D1<base_t>: public IALegendreIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         D1();

         /**
          * @brief Destructor
          */
         virtual ~D1() = default;

      protected:
         /**
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const override;

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
         virtual void applyPostOperator(Matrix& rOut, const bool isEven) const override;

         /**
          * @brief Apply pre FFT operator for component wise openerations
          *
          * @param in   Input values
          * @param useReal Real vs Imag flag
          */
         void applyPreOperator(const MatrixZ& in, const bool useReal, const bool isEven) const final;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          * @param useReal Real vs Imag flag
          */
         virtual void applyPostOperator(MatrixZ& rOut, const bool useReal, const bool isEven) const override;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_ALegendre_INTEGRATOR_BASE_D1_HPP
