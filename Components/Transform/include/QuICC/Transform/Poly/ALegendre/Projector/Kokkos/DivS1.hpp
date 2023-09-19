/**
 * @file DivS1.hpp
 * @brief Implementation of the associated Legendre based 1/sin P projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_DIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_DIVS1_HPP

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
   class DivS1;

   /**
    * @brief Implementation of the associated Legendre based P integrator
    */
   template <>
   class DivS1<kokkos_t>: public KokkosIALegendreProjector
   {
      public:
        /**
         * @brief Constructor
         */
        DivS1() = default;

        /**
         * @brief Destructor
         */
        virtual ~DivS1() = default;

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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_DIVS1_HPP
