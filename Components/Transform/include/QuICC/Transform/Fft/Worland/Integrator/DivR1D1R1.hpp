/**
 * @file DivR1D1R1.hpp
 * @brief Implementation of the Worland based DivR1D1R1 integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_DIVR1D1R1_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_DIVR1D1R1_HPP

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
#include "QuICC/Transform/Fft/Worland/Integrator/IWorlandIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   /**
    * @brief Implementation of the Worland based DivR1D1R1 integrator
    */
   class DivR1D1R1: public IWorlandIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         DivR1D1R1();

         /**
          * @brief Destructor
          */
         virtual ~DivR1D1R1();

      protected:
         /**
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const override;

         /**
          * @brief Initialise FFT backend
          */
         virtual void computeWorlandExpansion(const bool isEven) const override;

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyPreOperator(const Matrix& in, const bool isEven) const override;

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
         virtual void applyPreOperator(const MatrixZ& in, const bool useReal, const bool isEven) const override;

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

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_DIVR1D1R1_HPP
