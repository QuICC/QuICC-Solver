/**
 * @file NoSlip.hpp
 * @brief Implementation of the Bessel with NoSlip BC operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_NOSLIP_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_NOSLIP_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Polynomial/Bessel/NoSlip.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Integrator {

   /**
    * @brief Implementation of the Bessel with NoSlip BC operator
    */
   template <typename TOp, unsigned int Nbot = 0, unsigned int Ntop = 0>
   class NoSlip: public TOp
   {
      public:
         /**
          * @brief Constructor
          */
         NoSlip() = default;

         /**
          * @brief Destructor
          */
         virtual ~NoSlip() = default;

      protected:

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const override;
   };

   template <typename TOp,unsigned int Nbot, unsigned int Ntop> void NoSlip<TOp,Nbot,Ntop>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      this->template makeOperatorImpl<Polynomial::Bessel::NoSlip>(op, igrid, iweights, i);

      if constexpr(Nbot > 0)
      {
         op.rightCols(Nbot).setZero();
      }

      if constexpr(Ntop > 0)
      {
         op.leftCols(Ntop).setZero();
      }
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_NOSLIP_HPP
