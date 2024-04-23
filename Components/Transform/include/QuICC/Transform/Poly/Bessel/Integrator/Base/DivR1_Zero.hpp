/**
 * @file DivR1_Zero.hpp
 * @brief Implementation of the Bessel based 1/R1 integrator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_DIVR1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_DIVR1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/IBesselIntegrator.hpp"
#include "QuICC/Polynomial/Bessel/r_1SphJnl.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Integrator {

   template <class Impl>
   class DivR1_Zero;

   /**
    * @brief Implementation of the Bessel based 1/R1 integrator
    */
   template <>
   class DivR1_Zero<base_t>: public IBesselIntegrator
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
         TBc<Polynomial::Bessel::r_1SphJnl> r_1Jnl;
         r_1Jnl.template compute<Internal::MHDFloat>(op, nPoly, l, igrid, iweights);

         assert(op.rows() == igrid.size());
         assert(op.cols() == nPoly);
      }
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_INTEGRATOR_BASE_DIVR1_ZERO_HPP
