/**
 * @file Ll.hpp
 * @brief Implementation of the associated Legendre based l(l+1) P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_LL_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_LL_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/P.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class Ll;


   /**
    * @brief Implementation of the associated Legendre based l(l+1) P integrator
    */
   template <>
   class Ll<base_t>: public P<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         Ll() = default;

         /**
          * @brief Destructor
          */
         virtual ~Ll() = default;

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
          * @brief Storage for l(l+1) factors
          */
         mutable Array mLl;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_LL_HPP
