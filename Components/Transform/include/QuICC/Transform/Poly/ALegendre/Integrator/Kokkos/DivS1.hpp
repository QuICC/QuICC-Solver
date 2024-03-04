/**
 * @file PDivS1.hpp
 * @brief Implementation of the parallel associated Legendre based 1/Sin integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_DIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_DIVS1_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Tags.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/KokkosIALegendreIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class DivS1;

   /**
    * @brief Implementation of the associated Legendre based P integrator
    */
   template <>
   class DivS1<kokkos_t>: public KokkosIALegendreIntegrator
   {
      public:
         virtual void applyUnitOperator(const OpMatrixLZ &rOut,
           const OpMatrixLZ &in, const OpVectorI &scan,
           const int totalOpsCols) const override;

      protected:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const override;

   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_DIVS1_HPP
