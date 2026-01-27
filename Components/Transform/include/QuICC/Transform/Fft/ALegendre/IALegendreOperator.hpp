/**
 * @file IALegendreOperator.hpp
 * @brief Interface for a generic ALegendre FFT based operator
 */

#ifndef QUICC_TRANSFORM_FFT_ALEGENDRE_IALEGENDREOPERATOR_HPP
#define QUICC_TRANSFORM_FFT_ALEGENDRE_IALEGENDREOPERATOR_HPP

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
#include "QuICC/Transform/Fft/ALegendre/Setup.hpp"
#include "QuICC/Transform/Fft/IFftOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace ALegendre {

   /**
    * @brief Interface for a generic ALegendre FFT based operator
    */
   class IALegendreOperator: public IFftOperator
   {
      public:
         /// Typedef for the configuration class
         typedef Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         IALegendreOperator();

         /**
          * @brief Destructor
          */
         virtual ~IALegendreOperator();

         /**
          * @brief Initialise the transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedTransformSetup spSetup) const override;

         /**
          * @brief Initialise the transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const override;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const Matrix& in) const override;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const MatrixZ& in) const override;

      protected:
         /**
          * @brief Polynomial setup object providing the sizes
          */
         mutable SharedSetup    mspSetup;

         /**
          * @brief Compute transform R2C (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const Matrix& in) const override;
   };

}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_IWORLANDOPERATOR_HPP
