/**
 * @file LlD1.hpp
 * @brief Parallel Implementation of the associated Legendre based l(l+1) D projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_LlD1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_LlD1_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/D1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template <class Impl>
   class LlD1;

   /**
    * @brief Implementation of the associated Legendre based l(l+1) D projector
    */
   template <>
   class LlD1<kokkos_t>: public D1<kokkos_t>
   {
      public:
        /**
         * @brief Constructor
         */
        LlD1() = default;

        /**
         * @brief Destructor
         */
        virtual ~LlD1() = default;

      private:
        virtual void applyUnitOperator(const OpMatrixLZ &rOut,
           const OpMatrixLZ &in, const OpVectorI &scan,
           const int totalOpsCols) const;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const;

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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_LlD1_HPP
