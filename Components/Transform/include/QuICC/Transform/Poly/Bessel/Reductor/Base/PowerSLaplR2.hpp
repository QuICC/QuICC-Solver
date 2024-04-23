/**
 * @file PowerSLaplR2.hpp
 * @brief Implementation of the Bessel based spherical Laplacian R^2 power spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_POWERSLAPLR2_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_POWERSLAPLR2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/IBesselPower.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"
#include "QuICC/Polynomial/Bessel/slaplSphJnl.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   template <class Impl>
   class PowerSLaplR2;

   /**
    * @brief Implementation of the Bessel based Spherical Laplacian R^2 power spectrum operator
    */
   template <>
   class PowerSLaplR2<base_t>: public Bessel::Reductor::IBesselPower
   {
      public:
         /**
          * @brief Constructor
          */
         PowerSLaplR2();

         /**
          * @brief Destructor
          */
         virtual ~PowerSLaplR2() = default;

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const override;

         /**
          * @brief Make operator
          */
         template <template <typename> class TBc> void makeOperatorImpl(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const;
      private:
   };

   template <template <typename> class TBc> void PowerSLaplR2<base_t>::makeOperatorImpl(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fastSize(i);

      // Build operator
      op.resize(igrid.size(), nPoly);
      TBc<Polynomial::Bessel::slaplSphJnl> bjnl;
      bjnl.template compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array());

      TBc<Polynomial::Bessel::SphJnl> fjnl;

      eop.resize(igrid.size(), nPoly);
      fjnl.template compute<MHDFloat>(eop, nPoly, l, igrid, iweights);
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_POWERSLAPLR2_HPP
