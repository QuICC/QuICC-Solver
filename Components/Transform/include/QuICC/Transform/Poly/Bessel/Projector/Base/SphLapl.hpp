/**
 * @file SphLapl.hpp
 * @brief Implementation of the Bessel based spherical laplacian projector
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_SPHLAPL_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_SPHLAPL_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/IBesselProjector.hpp"
#include "QuICC/Polynomial/Bessel/slaplSphJnl.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Projector {

   template <class Impl>
   class SphLapl;

   /**
    * @brief Implementation of the Bessel based spherical laplacian projector
    */
   template <>
   class SphLapl<base_t>: public Bessel::Projector::IBesselProjector
   {
      public:
         /**
          * @brief Constructor
          */
         SphLapl();

         /**
          * @brief Destructor
          */
         ~SphLapl() = default;

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

   template <template <typename> class TBc> void SphLapl<base_t>::makeOperatorImpl(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      Internal::Matrix top(igrid.size(), nPoly);
      TBc<Polynomial::Bessel::slaplSphJnl> jnl;
      jnl.template compute<Internal::MHDFloat>(top, nPoly, l, igrid, Internal::Array());

      op = top.cast<MHDFloat>();
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_SPHLAPL_HPP
