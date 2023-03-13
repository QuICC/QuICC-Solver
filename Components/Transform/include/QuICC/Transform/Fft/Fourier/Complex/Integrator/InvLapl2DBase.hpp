/**
 * @file InvLapl2D.hpp
 * @brief Implementation of the Fourier based inverse 2D Laplacian integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_INVLAPL2DBASE_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_INVLAPL2DBASE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/IComplexIntegrator.hpp"
#include "QuICC/Transform/Fft/Fourier/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {
   template<class>
   class InvLapl2D;

   /**
    * @brief Implementation of the Fourier based inverse 2D Laplacian integrator
    */
   template<>
   class InvLapl2D<base_t>: public IComplexIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         InvLapl2D() = default;

         /**
          * @brief Destructor
          */
         ~InvLapl2D() = default;

      protected:

      private:
         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         void applyPostOperator(MatrixZ& rOut) const final;
   };

} // namespace Integrator
} // namespace Complex
} // namespace Fourier
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_INVLAPL2DBASE_HPP
