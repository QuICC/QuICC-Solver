/**
 * @file DivR1_Zero.hpp
 * @brief Implementation of the Bessel based 1/R projector but 0 mode is zeroed
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_DIVR1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_DIVR1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/IBesselProjector.hpp"
#include "QuICC/Polynomial/Bessel/r_1SphJnl.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Projector {

   template <class Impl>
   class DivR1_Zero;

   /**
    * @brief Implementation of the Bessel based 1/R projector but 0 mode is zeroed
    */
   template <>
   class DivR1_Zero<base_t>: public Bessel::Projector::IBesselProjector
   {
      public:
         /**
          * @brief Constructor
          */
         DivR1_Zero();

         /**
          * @brief Destructor
          */
         ~DivR1_Zero() = default;

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

   template <template <typename> class TBc> void DivR1_Zero<base_t>::makeOperatorImpl(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
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
         TBc<Polynomial::Bessel::r_1SphJnl> jnl;
         jnl.template compute<Internal::MHDFloat>(top, nPoly, l, igrid, Internal::Array());

         op = top.cast<MHDFloat>();
      }
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_DIVR1_ZERO_HPP
