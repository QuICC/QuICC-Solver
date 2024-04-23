/**
 * @file R1_Zero.hpp
 * @brief Implementation of the Bessel based R1_Zero integrator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_R1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_R1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/IBesselIntegrator.hpp"
#include "QuICC/Polynomial/Bessel/rSphJnl.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Integrator {

   template <class Impl>
   class R1_Zero;

   /**
    * @brief Implementation of the Bessel based R1_Zero integrator
    */
   template <>
   class R1_Zero<base_t>: public IBesselIntegrator
   {
      public:
         /**
          * @brief Constructor
          */
         R1_Zero();

         /**
          * @brief Destructor
          */
         ~R1_Zero() = default;

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

   template <template <typename> class TBc> void R1_Zero<base_t>::makeOperatorImpl(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
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
         TBc<Polynomial::Bessel::rSphJnl> jnl;
         jnl.template compute<Internal::MHDFloat>(op, nPoly, l, igrid, iweights);
      }
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_R1_ZERO_HPP
