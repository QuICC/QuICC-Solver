/**
 * @file IReductor.hpp
 * @brief Interface for a Worland based reduction operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_IREDUCTOR_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_IREDUCTOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/IOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   /**
    * @brief Interface for a Worland based energy operator
    */
   class IReductor: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IReductor();

         /**
          * @brief Destructor
          */
         virtual ~IReductor();

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_IREDUCTOR_HPP
