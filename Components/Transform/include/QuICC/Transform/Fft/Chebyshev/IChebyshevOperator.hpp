/**
 * @file IChebyshevOperator.hpp
 * @brief Interface for a generic Chebyshev FFT based operator
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_ICHEBYSHEVOPERATOR_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_ICHEBYSHEVOPERATOR_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"
#include "QuICC/Transform/Fft/IFftOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

   /**
    * @brief Interface for a generic Chebyshev FFT based operator
    */
   class IChebyshevOperator: public IFftOperator
   {
      public:
         /// Typedef for the configuration class
         typedef Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         IChebyshevOperator();

         /**
          * @brief Destructor
          */
         virtual ~IChebyshevOperator();

         /**
          * @brief Initialise the transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedSetupType spSetup) const;

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

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_ICHEBYSHEVOPERATOR_HPP
