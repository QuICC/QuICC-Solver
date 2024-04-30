/**
 * @file D1.hpp
 * @brief Implementation of the Bessel based D1 projector
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_D1_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_D1_HPP

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Polynomial/Bessel/dSphJnl.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/IBesselProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Projector {

   template <class Impl>
   class D1;

   /**
    * @brief Implementation of the Bessel based P projector
    */
   template <>
   class D1<base_t>: public Bessel::Projector::IBesselProjector
   {
      public:
         /**
          * @brief Constructor
          */
         D1();

         /**
          * @brief Destructor
          */
         ~D1() = default;

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

   template <template <typename> class TBc> void D1<base_t>::makeOperatorImpl(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      Internal::Matrix top(igrid.size(), nPoly);
      TBc<Polynomial::Bessel::dSphJnl> jnl;
      jnl.template compute<Internal::MHDFloat>(top, nPoly, l, igrid, Internal::Array());

      op = top.cast<MHDFloat>();
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_D1_HPP
