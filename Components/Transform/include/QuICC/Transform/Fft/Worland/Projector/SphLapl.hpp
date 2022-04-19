/** 
 * @file SphLapl.hpp
 * @brief Implementation of the Worland based spherial laplacian projector
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_SPHLAPL_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_SPHLAPL_HPP

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
    * @brief Implementation of the Worland based P projector
    */ 
   class SphLapl: public IWorlandProjector
   {
      public:
         /**
          * @brief Constructor
          */
         SphLapl();

         /**
          * @brief Destructor
          */
         virtual ~SphLapl();
         
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

#endif // QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_SPHLAPL_HPP
