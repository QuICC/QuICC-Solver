/**
 * @file P_Zero.hpp
 * @brief Implementation of the Bessel based integrator, zero for l = 0
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_P_ZERO_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_P_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/IBesselIntegrator.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Integrator {

   template <class Impl>
   class P_Zero;

   /**
    * @brief Implementation of the Bessel based integrator, zero for l = 0
    */
   template <>
   class P_Zero<base_t>: public IBesselIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         P_Zero();

         /**
          * @brief Destructor
          */
         ~P_Zero() = default;

      protected:
         /**
          * @brief Make operator
          */
         template <template <typename> class TBc> void makeOperatorImpl(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const;

      private:
         /**
          * @brief Apply ith operator
          */
         void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const final;
   };

   template <template <typename> class TBc> void P_Zero<base_t>::makeOperatorImpl(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      if(l == 0)
      {
         op.setZero();
      }
      else
      {
         Internal::Matrix top(igrid.size(), nPoly);
         TBc<Polynomial::Bessel::SphJnl> jnl;
         jnl.template compute<Internal::MHDFloat>(top, nPoly, l, igrid, iweights);

         op = top.cast<MHDFloat>();
      }
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_P_ZERO_HPP
