/**
 * @file DivS1Dp.hpp
 * @brief Implementation of the associated Legendre based 1/sin P d_phi projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_DIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_DIVS1DP_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/DivS1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template <class Impl>
   class DivS1Dp;

   /**
    * @brief Implementation of the associated Legendre based 1/sin P d_phi projector
    */
   template <>
   class DivS1Dp<kokkos_t>: public DivS1<kokkos_t>
   {
      public:
        /**
         * @brief Constructor
         */
        DivS1Dp() = default;

        /**
         * @brief Destructor
         */
        virtual ~DivS1Dp() = default;

      private:
        virtual void applyUnitOperator(const OpMatrixLZ &rOut,
           const OpMatrixLZ &in, const OpVectorI &scan,
           const int totalOpsCols) const;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_DIVS1DP_HPP
