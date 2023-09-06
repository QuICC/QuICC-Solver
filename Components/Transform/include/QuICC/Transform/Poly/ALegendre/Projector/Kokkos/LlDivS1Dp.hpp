/**
 * @file LlDivS1Dp.hpp
 * @brief Parallel Implementation of the associated Legendre based 1/sin l(l+1) P d_phi P projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_LLDIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_LLDIVS1DP_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
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
   class LlDivS1Dp;

   /**
    * @brief Implementation of the associated Legendre based 1/sin l(l+1) P d_phi projector
    */
   template <>
   class LlDivS1Dp<kokkos_t>: public DivS1<kokkos_t>
   {
      public:/**
         * @brief Constructor
         */
        LlDivS1Dp() = default;

        /**
         * @brief Destructor
         */
        virtual ~LlDivS1Dp() = default;

      private:
        virtual void applyUnitOperator(const OpMatrixLZ &rOut,
           const OpMatrixLZ &in, const OpVectorI &scan,
           const int totalOpsCols) const override;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const override;

        /**
         * @brief l(l+1) scaling factors
         */
         virtual void initSpecial() const;

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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_LLDIVS1DP_HPP
