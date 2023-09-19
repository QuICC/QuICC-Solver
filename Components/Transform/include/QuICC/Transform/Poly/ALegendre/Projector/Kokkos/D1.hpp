/**
 * @file D1.hpp
 * @brief Implementation of the associated Legendre based D projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_D1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_D1_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Tags.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/KokkosIALegendreProjector.hpp"


namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template <class Impl>
   class D1;

   /**
    * @brief Implementation of the associated Legendre based D projector
    */
   template <>
   class D1<kokkos_t>: public KokkosIALegendreProjector
   {
      public:
        /**
         * @brief Constructor
         */
        D1() = default;

        /**
         * @brief Destructor
         */
        virtual ~D1() = default;

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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_D1_HPP
