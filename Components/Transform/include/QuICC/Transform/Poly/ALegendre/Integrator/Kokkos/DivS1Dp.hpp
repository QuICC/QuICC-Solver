/**
 * @file DivS1Dp.hpp
 * @brief Implementation of the associated Legendre based 1/Sin D_phi parallel integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_DIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_DIVS1DP_HPP

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
   class DivS1Dp;

   /**
    * @brief Implementation of the associated Legendre based P integrator
    */
   template <>
   class DivS1Dp<kokkos_t>: public DivS1<kokkos_t>
   {
      private:
        virtual void applyUnitOperator(const OpMatrixLZ &rOut,
           const OpMatrixLZ &in, const OpVectorI &scan,
           const int totalOpsCols) const;

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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_DIVS1DP_HPP
