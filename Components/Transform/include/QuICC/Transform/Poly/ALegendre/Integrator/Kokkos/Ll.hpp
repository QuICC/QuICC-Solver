/**
 * @file Ll.hpp
 * @brief Implementation of the associated Legendre based LL integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_Ll_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_Ll_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/P.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class Ll;

   /**
    * @brief Implementation of the associated Legendre based P integrator
    */
   template <>
   class Ll<kokkos_t>: public P<kokkos_t>
   {
      private:
        virtual void applyUnitOperator(const OpMatrixLZ &rOut,
           const OpMatrixLZ &in, const OpVectorI &scan,
           const int totalOpsCols) const;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const;

        /**
         * @brief l(l+1) scaling factors
         */
        virtual void initSpecial() const override;

        /**
         * @brief Storage for l(l+1) factors
         */
        mutable Array mLl;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_Ll_HPP
