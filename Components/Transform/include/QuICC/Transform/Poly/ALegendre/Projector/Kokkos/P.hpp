/**
 * @file PP.hpp
 * @brief Implementation of the associated Legendre based PP Projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_PP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_PP_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Tags.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/KokkosIALegendreProjector.hpp"



namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template <class Impl>
   class P;

   /**
    * @brief Implementation of the associated Legendre based PP Projector
    */
   template <>
   class P<kokkos_t>: public KokkosIALegendreProjector
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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_P_HPP
