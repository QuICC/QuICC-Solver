/**
 * @file LlDivS1Dp.hpp
 * @brief Implementation of the associated Legendre based l(l+1)/Sin D_phi integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_LLDIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_LLDIVS1DP_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/DivS1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class LlDivS1Dp;

   /**
    * @brief Implementation of the associated Legendre based P integrator
    */
   template <>
   class LlDivS1Dp<kokkos_t>: public DivS1<kokkos_t>
   {
      private:
        virtual void applyUnitOperator(const OpMatrixLZ &rOut,
           const OpMatrixLZ &in, const OpVectorI &scan,
           const int totalOpsCols) const override;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const override;

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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_LLDIVS1DP_HPP
