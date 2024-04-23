/**
 * @file PowerR2.hpp
 * @brief Implementation of the Bessel based R^2 power spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_POWERR2_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_POWERR2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/IBesselPower.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   template <class Impl>
   class PowerR2;

   /**
    * @brief Implementation of the Bessel based R^2 power spectrum operator
    */
   template <>
   class PowerR2<base_t>: public Bessel::Reductor::IBesselPower
   {
      public:
         /**
          * @brief Constructor
          */
         PowerR2();

         /**
          * @brief Destructor
          */
         virtual ~PowerR2() = default;

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const override;

         /**
          * @brief Make operator
          */
         template <template <typename> class TBc > void makeOperatorImpl(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const;
      private:
   };

   template <template <typename> class TBc> void PowerR2<base_t>::makeOperatorImpl(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fastSize(i);

      // Build operator
      op.resize(igrid.size(), nPoly);
      TBc<Polynomial::Bessel::SphJnl> bjnl;
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

#endif // QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_POWERR2_HPP
