/**
 * @file P.hpp
 * @brief Implementation of the Bessel based P projector
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_P_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_P_HPP

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/IBesselProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Projector {

   template <class Impl>
   class P;

   /**
    * @brief Implementation of the Bessel based P projector
    */
   template <>
   class P<base_t>: public Bessel::Projector::IBesselProjector
   {
      public:
         /**
          * @brief Constructor
          */
         P();

         /**
          * @brief Destructor
          */
         ~P() = default;

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

   template <template <typename> class TBc> void P<base_t>::makeOperatorImpl(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);

      // Build operator
      int nPoly = this->mspSetup->fastSize(i);
      op.resize(igrid.size(), nPoly);
      TBc<Polynomial::Bessel::SphJnl> jnl;
      jnl.template compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array());
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_P_HPP
