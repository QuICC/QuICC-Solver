/** 
 * @file DivR1_Zero.hpp
 * @brief Implementation of the Worland based 1/R projector and zero for l = 0
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_DIVR1_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_DIVR1_ZERO_HPP

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
#include "QuICC/Transform/Fft/Worland/Projector/IWorlandProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   /**
    * @brief Implementation of the Worland based 1/R projector and zero for l = 0
    */ 
   class DivR1_Zero: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         DivR1_Zero();

         /**
          * @brief Destructor
          */
         virtual ~DivR1_Zero();
         
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
         virtual void applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const override;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          * @param useReal Real vs Imag flag
          */
         virtual void applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const override;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_DIVR1_ZERO_HPP
