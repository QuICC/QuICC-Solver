/**
 * @file P.hpp
 * @brief Implementation of the parallel associated Legendre based D integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_D1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_D1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Tags.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/KokkosIALegendreIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class D1;

   /**
    * @brief Implementation of the associated Legendre based P integrator
    */
   template <>
   class D1<kokkos_t>: public KokkosIALegendreIntegrator
   {
      public:

        virtual void applyUnitOperator(const OpMatrixLZ &rOut,
           const OpMatrixLZ &in, const OpVectorI &scan,
           const int totalOpsCols) const;

      protected:
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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_D1_HPP
