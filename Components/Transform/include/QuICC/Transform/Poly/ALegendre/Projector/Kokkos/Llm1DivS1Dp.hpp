/**
 * @file Llm1DivS1Dp.hpp
 * @brief Parallel Implementation of the associated Legendre based 1/sin [l(l+1)]-1 P d_phi P projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_LLM1DIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_LLM1DIVS1DP_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/DivS1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template <class Impl>
   class Llm1DivS1Dp;

   /**
    * @brief Implementation of the associated Legendre based 1/sin [l(l+1)-1] P d_phi projector
    */
   template <>
   class Llm1DivS1Dp<kokkos_t>: public DivS1<kokkos_t>
   {
      public:/**
         * @brief Constructor
         */
        Llm1DivS1Dp() = default;

        /**
         * @brief Destructor
         */
        virtual ~Llm1DivS1Dp() = default;

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
         virtual void initSpecial() const;

         /**
          * @brief Storage for l(l+1)-1 factors
          */
         mutable Array mLlm1;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_LLM1DIVS1DP_HPP
