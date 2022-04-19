/**
 * @file IChebyshevReductor.hpp
 * @brief Interface for a generic Chebyshev FFT based reductor, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_ICHEBYSHEVREDUCTOR_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_ICHEBYSHEVREDUCTOR_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

   /**
    * @brief Interface for a generic Chebyshev FFT based reductor
    */
   class IChebyshevReductor: public IChebyshevOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IChebyshevReductor();

         /**
          * @brief Destructor
          */
         virtual ~IChebyshevReductor();

      protected:

      private:
   };

}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_ICHEBYSHEVREDUCTOR_HPP
