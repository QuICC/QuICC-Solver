/**
 * @file Insulating.hpp
 * @brief Implementation of the Bessel with Insulating BC operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_INSULATING_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_INSULATING_HPP

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Polynomial/Bessel/Insulating.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Integrator {

   /**
    * @brief Implementation of the Bessel with Insulating BC operator
    */
   template <typename TOp, unsigned int N = 0>
   class Insulating: public TOp
   {
      public:
         /**
          * @brief Constructor
          */
         Insulating() = default;

         /**
          * @brief Destructor
          */
         virtual ~Insulating() = default;

      protected:

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const override;
   };

   template <typename TOp, unsigned int N> void Insulating<TOp,N>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      this->template makeOperatorImpl<Polynomial::Bessel::Insulating>(op, igrid, iweights, i);

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

#endif // QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_INSULATING_HPP
