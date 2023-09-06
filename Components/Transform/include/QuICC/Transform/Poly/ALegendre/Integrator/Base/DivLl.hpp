/**
 * @file DivLl.hpp
 * @brief Implementation of the associated Legendre based 1/l(l+1) P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_DIVLL_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_DIVLL_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/P.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class DivLl;

   /**
    * @brief Implementation of the associated Legendre based 1/(l(l+1)) P integrator
    */
   template <>
   class DivLl<base_t>: public P<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         DivLl() = default;

         /**
          * @brief Destructor
          */
         virtual ~DivLl() = default;

      protected:

      private:
         /**
          * @brief Apply ith operator
          */
         void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const final;

         /**
          * @brief l(l+1) scaling factors
          */
         void initSpecial() const;

         /**
          * @brief Storage for 1/l(l+1) factors
          */
         mutable Array mDivLl;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_DIVLL_HPP
