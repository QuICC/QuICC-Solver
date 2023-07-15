/**
 * @file DivLlDivS1.hpp
 * @brief Implementation of the associated Legendre based 1/l(l+1) 1/Sin parallel integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_DIVLLDIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_DIVLLDIVS1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/DivS1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class DivLlDivS1;

   /**
    * @brief Implementation of the associated Legendre based P integrator
    */
   template <>
   class DivLlDivS1<kokkos_t>: public DivS1<kokkos_t>
   {
      private:
        virtual void applyUnitOperator(const OpMatrixLZ &rOut,
           const OpMatrixLZ &in, const OpVectorI &scan,
           const int totalOpsCols) const;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const;

        /**
         * @brief l(l+1) scaling factors
         */
        virtual void initSpecial() const override;

        /**
         * @brief Storage for l(l+1) factors
         */
        mutable Array mDivLl;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_DIVLLDIVS1_HPP
