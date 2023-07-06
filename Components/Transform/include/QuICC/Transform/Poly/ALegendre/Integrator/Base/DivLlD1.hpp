/**
 * @file DivLlD1.hpp
 * @brief Implementation of the associated Legendre based 1/l(l+1) D integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_DIVLLD1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_DIVLLD1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/D1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template <class Impl>
   class DivLlD1;

   /**
    * @brief Implementation of the associated Legendre based 1/(l(l+1)) D integrator
    */
   template <>
   class DivLlD1<base_t>: public D1<base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         DivLlD1() = default;

         /**
          * @brief Destructor
          */
         virtual ~DivLlD1() = default;

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
          * @brief Storage for 1/l(l+1) factors
          */
         mutable Array mDivLl;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_BASE_DIVLLD1_HPP
