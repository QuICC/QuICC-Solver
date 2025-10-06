/**
 * @file DivY1D2Y1.hpp
 * @brief Implementation of the Chebyshev based radial part of spherical Laplacian 1/Y D2 Y D projector, with linear map y = ax + b
 * It is equivalent to SphRadLapl  
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_DIVY1D2Y1_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_DIVY1D2Y1_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/IChebyshevProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   /**
    * @brief Implementation of the Chebyshev based radial part of spherical Laplacian 1/Y^2 D Y^2 D projector, with linear map y = ax + b
    */
   class DivY1D2Y1: public IChebyshevProjector
   {
      public:
         /**
          * @brief Constructor
          */
         DivY1D2Y1();

         /**
          * @brief Destructor
          */
         virtual ~DivY1D2Y1();

      protected:
         /**
          * @brief Initialize storage
          */
         virtual void initBackend() const override;

      private:
         /**
          * @brief Initialize solver operators
          */
         virtual void initOperator() const override;

         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyPreOperator(Matrix& rOut, const Matrix& in) const override;

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
         virtual void applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const override;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          * @param useReal Real vs Imag flag
          */
         virtual void applyPostOperator(MatrixZ& rOut, const Matrix& tmp,  const bool useReal) const override;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_DIVY1D2Y1_HPP
