/**
 * @file IWorlandReductor.hpp
 * @brief Interface for a Worland based reduction operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_IWORLANDREDUCTOR_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_IWORLANDREDUCTOR_HPP

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
#include "QuICC/Transform/Poly/Worland/IWorlandOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   /**
    * @brief Interface for a Worland based energy operator
    */
   class IWorlandReductor: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IWorlandReductor();

         /**
          * @brief Destructor
          */
         virtual ~IWorlandReductor();

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
         mutable internal::Array  mGrid;

         /**
          * @brief Storage for the quadrature weights
          */
         mutable internal::Array  mWeights;

      private:
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_IWORLANDREDUCTOR_HPP
