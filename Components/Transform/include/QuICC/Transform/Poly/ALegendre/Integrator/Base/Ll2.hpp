/**
 * @file Ll2.hpp
 * @brief Implementation of the associated Legendre based l(l+1)^2 P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_LL2_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_LL2_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/P.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class Ll2;

   /**
    * @brief Implementation of the associated Legendre based l(l+1)^2 P integrator
    */
   template <>
   class Ll2<base_t>: public P<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         Ll2() = default;

         /**
          * @brief Destructor
          */
         virtual ~Ll2() = default;

      protected:

      private:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const override;

         /**
          * @brief l(l+1) scaling factors
          */
         virtual void initSpecial() const override;

         /**
          * @brief Storage for l(l+1)^2 factors
          */
         mutable Array mLl2;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_LL2_HPP
