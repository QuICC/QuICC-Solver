/**
 * @file DivR1D1R1_Zero.hpp
 * @brief Implementation of the Bessel based 1/R D R projector and zero for l = 0
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_DIVR1D1R1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_DIVR1D1R1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/IBesselProjector.hpp"
#include "QuICC/Polynomial/Bessel/r_1drSphJnl.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Projector {

   template <class Impl>
   class DivR1D1R1_Zero;

   /**
    * @brief Implementation of the Bessel based 1/R D R projector and zero for l = 0
    */
   template <>
   class DivR1D1R1_Zero<base_t>: public Bessel::Projector::IBesselProjector
   {
      public:
         /**
          * @brief Constructor
          */
         DivR1D1R1_Zero();

         /**
          * @brief Destructor
          */
         ~DivR1D1R1_Zero() = default;

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

   template <template <typename> class TBc> void DivR1D1R1_Zero<base_t>::makeOperatorImpl(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
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
         TBc<Polynomial::Bessel::r_1drSphJnl<QuICC::Polynomial::Bessel::recurrence_t>> jnl;
         jnl.template compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array());
      }
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_PROJECTOR_BASE_DIVR1D1R1_ZERO_HPP
