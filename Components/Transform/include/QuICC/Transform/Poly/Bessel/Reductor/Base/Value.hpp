/**
 * @file Value.hpp
 * @brief Implementation of the Bessel with Value BC operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_VALUE_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_VALUE_HPP

// External includes
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

namespace Reductor {

   /**
    * @brief Implementation of the Bessel with Value BC operator
    */
   template <typename TOp>
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
         virtual void makeOperator(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const override;
   };

   template <typename TOp> void Value<TOp>::makeOperator(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      this->template makeOperatorImpl<Polynomial::Bessel::Value>(op, eop, igrid, iweights, i);
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_VALUE_HPP
