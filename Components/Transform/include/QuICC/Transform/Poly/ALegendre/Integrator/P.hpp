/**
 * @file P.hpp
 * @brief Implementation of the associated Legendre based P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_P_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_P_HPP

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
#include "QuICC/Transform/Poly/ALegendre/Integrator/IALegendreIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   /**
    * @brief Implementation of the associated Legendre based P integrator
    */
   template<typename OpTypes = IALegendreOperatorTypes>
   class P: public IALegendreIntegrator<OpTypes>
   {
      public:
        // using Array = typename OpTypes::Array;
        using OpArray = typename OpTypes::OpArray;
        using OpMatrix = typename OpTypes::OpMatrix;
        using OpMatrixR = typename OpTypes::OpMatrixR;
        using OpMatrixCR = typename OpTypes::OpMatrixCR;

        /**
         * @brief Constructor
         */
        P();

        /**
         * @brief Destructor
         */
        virtual ~P();

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const;

       private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const;

   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_P_HPP
