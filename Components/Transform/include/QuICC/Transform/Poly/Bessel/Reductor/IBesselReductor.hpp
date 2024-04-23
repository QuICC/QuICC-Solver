/**
 * @file IBesselReductor.hpp
 * @brief Interface for a Bessel based reduction operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_IBESSELREDUCTOR_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_IBESSELREDUCTOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/Bessel/IBesselOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   /**
    * @brief Interface for a Bessel based energy operator
    */
   class IBesselReductor: public IBesselOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IBesselReductor();

         /**
          * @brief Destructor
          */
         virtual ~IBesselReductor() = default;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const;

      protected:
         /**
          * @brief Storage for the operators
          */
         mutable std::vector<Matrix>  mOps;

         /**
          * @brief Storage for the quadrature grid
          */
         mutable Internal::Array  mGrid;

         /**
          * @brief Storage for the quadrature weights
          */
         mutable Internal::Array  mWeights;

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_IBESSELREDUCTOR_HPP
