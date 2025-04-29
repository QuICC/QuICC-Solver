/**
 * @file P.hpp
 * @brief Implementation of the associated Legendre based P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_P_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_P_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Tags.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/IALegendreIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

    template <class Impl>
    class P;

   /**
    * @brief Implementation of the associated Legendre based P integrator
    */
   template <>
   class P<base_t>: public IALegendreIntegrator
   {
      public:

        /**
         * @brief Constructor
         */
        P() = default;

        /**
         * @brief Destructor
         */
        virtual ~P() = default;

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const;

       private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const;

   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_P_HPP
