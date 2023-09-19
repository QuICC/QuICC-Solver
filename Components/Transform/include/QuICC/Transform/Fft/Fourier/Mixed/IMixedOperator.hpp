/**
 * @file IMixedOperator.hpp
 * @brief Interface for a generic Mixed FFT based operator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_IMIXEDOPERATOR_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_IMIXEDOPERATOR_HPP

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
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"
#include "QuICC/Transform/Fft/IFftOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

   /**
    * @brief Interface for a generic Mixed FFT based operator
    */
   class IMixedOperator: public IFftOperator
   {
      public:
         /// Typedef for the configuration class
         typedef Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         IMixedOperator();

         /**
          * @brief Destructor
          */
         virtual ~IMixedOperator();

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
         void init(SharedTransformSetup spSetup, const internal::Array& igrid, const internal::Array& iweights) const override;

      protected:
         /**
          * @brief Polynomial setup object providing the sizes
          */
         mutable SharedSetup    mspSetup;

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_IMIXEDOPERATOR_HPP
