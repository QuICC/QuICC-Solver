/**
 * @file Value.hpp
 * @brief Implementation of the Bessel with Value BC operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_VALUE_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_VALUE_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Polynomial/Bessel/Value.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Integrator {

   /**
    * @brief Implementation of the Bessel with Value BC operator
    */
   template <typename TOp, unsigned int N = 0>
   class Value: public TOp
   {
      public:
         /**
          * @brief Constructor
          */
         Value() = default;

         /**
          * @brief Destructor
          */
         virtual ~Value() = default;

      protected:

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const override;
   };

   template <typename TOp,unsigned int N> void Value<TOp,N>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      this->template makeOperatorImpl<Polynomial::Bessel::Value>(op, igrid, iweights, i);

      if constexpr(N > 0)
      {
         op.bottomRows(N).setZero();
      }
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_VALUE_HPP
