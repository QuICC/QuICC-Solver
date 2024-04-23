/**
 * @file PowerD1R1.hpp
 * @brief Implementation of the Bessel based D R power spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_POWERD1R1_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_POWERD1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Tags.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/IBesselPower.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"
#include "QuICC/Polynomial/Bessel/r_1drSphJnl.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   template <class Impl>
   class PowerD1R1;

   /**
    * @brief Implementation of the Bessel based D R power spectrum operator
    */
   template <>
   class PowerD1R1<base_t>: public Bessel::Reductor::IBesselPower
   {
      public:
         /**
          * @brief Constructor
          */
         PowerD1R1();

         /**
          * @brief Destructor
          */
         virtual ~PowerD1R1() = default;

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

   template <template <typename> class TBc> void PowerD1R1<base_t>::makeOperatorImpl(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int l = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fastSize(i);

      // Build operator
      op.resize(igrid.size(), nPoly);
      TBc<Polynomial::Bessel::r_1drSphJnl<QuICC::Polynomial::Bessel::recurrence_t>> bjnl;
      bjnl.template compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array());

      eop.resize(igrid.size(), nPoly);
      TBc<Polynomial::Bessel::SphJnl> fjnl;
      if(l == 0)
      {
         eop.setZero();
      }
      else
      {
         fjnl.template compute<MHDFloat>(eop, nPoly, std::abs(l-1), igrid, iweights);
      }
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_BASE_POWERD1R1_HPP
