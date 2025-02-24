/**
 * @file I3QpDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I3QpDiags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I3QpDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   using namespace Internal::Literals;

   I3QpDiags::I3QpDiags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::I3QpDiags(alpha, MHD_MP(-0.5), l, q), mI2(alpha, l, 0)
   {
      if(q > 1)
      {
         throw std::logic_error("I3Qp: Truncation for q>1 is not implemented");
      }
   }

   I3QpDiags::ACoeff_t I3QpDiags::d_3(const ACoeff_t& n) const
   {
      ACoeff_t val;
      val = 8.0_mp*(l<1>() + n - 2.0_mp)*(l<1>() + n - 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n - 3.0_mp)/((l<1>() + 2.0_mp*n - 5.0_mp)*(l<1>() + 2.0_mp*n - 4.0_mp)*(l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp));

      return this->normalizeDiag(n,-3,1)*val;
   }

   I3QpDiags::ACoeff_t I3QpDiags::d_2(const ACoeff_t& n) const
   {
      ACoeff_t val;
      val = -4.0_mp*(l<1>() + n - 1.0_mp)*(8.0_mp*l<2>() + 4.0_mp*l<1>()*n - 2.0_mp*l<1>() - 4.0_mp*n*n + 8.0_mp*n + 5.0_mp)/((l<1>() + 2.0_mp*n - 4.0_mp)*(l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp));

      return this->normalizeDiag(n,-2,1)*val;
   }

   I3QpDiags::ACoeff_t I3QpDiags::d_1(const ACoeff_t& n) const
   {
      ACoeff_t val;
      val = 4.0_mp*(4.0_mp*l<3>() - 12.0_mp*l<2>()*n + 6.0_mp*l<2>() - 24.0_mp*l<1>()*n*n + 12.0_mp*l<1>()*n + 20.0_mp*l<1>() - 8.0_mp*n*n*n + 4.0_mp*n*n + 10.0_mp*n + 3.0_mp)/((l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp));

      // Correct if truncation q == 1
      this->correctQ1(val, n, -1);

      return this->normalizeDiag(n,-1,1)*val;
   }

   I3QpDiags::ACoeff_t I3QpDiags::d0(const ACoeff_t& n) const
   {
      ACoeff_t val;
      val = 2.0_mp*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)*(12.0_mp*l<2>()*n + 2.0_mp*l<2>() + 4.0_mp*l<1>()*n - 10.0_mp*l<1>() - 8.0_mp*n*n*n - 4.0_mp*n*n + 10.0_mp*n - 3.0_mp)/((l<1>() + n)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp));

      // Correct if truncation q == 1
      this->correctQ1(val, n, 0);

      return this->normalizeDiag(n,0,1)*val;
   }

   I3QpDiags::ACoeff_t I3QpDiags::d1(const ACoeff_t& n) const
   {
      ACoeff_t val;
      val = (2.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 3.0_mp)*(12.0_mp*l<1>()*n + 10.0_mp*l<1>() + 4.0_mp*n*n + 8.0_mp*n - 5.0_mp)/(2.0_mp*(l<1>() + n)*(l<1>() + n + 1.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp)*(l<1>() + 2.0_mp*n + 4.0_mp));

      // Correct if truncation q == 1
      this->correctQ1(val, n, 1);

      return this->normalizeDiag(n,1,1)*val;
   }

   I3QpDiags::ACoeff_t I3QpDiags::d2(const ACoeff_t& n) const
   {
      ACoeff_t val;
      val = (2.0_mp*n + 1.0_mp)*(2.0_mp*n + 3.0_mp)*(2.0_mp*n + 3.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 3.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 5.0_mp)/(4.0_mp*(l<1>() + n)*(l<1>() + n + 1.0_mp)*(l<1>() + n + 2.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp)*(l<1>() + 2.0_mp*n + 4.0_mp)*(l<1>() + 2.0_mp*n + 5.0_mp));

      // Correct if truncation q == 1
      this->correctQ1(val, n, 2);

      return this->normalizeDiag(n,2,1)*val;
   }

   I3QpDiags::ACoeff_t I3QpDiags::d3(const ACoeff_t& n) const
   {
      ACoeff_t val = ACoeff_t::Zero(n.size());

      // Correct if truncation q == 1
      this->correctQ1(val, n, 3);

      return this->normalizeDiag(n,3,1)*val;
   }

   void I3QpDiags::correctQ1(ACoeff_t& val, const ACoeff_t& n, const int k) const
   {
      // Index where to apply correction in val
      auto i_ = val.size() - (k+2);

      // Only correct if truncation q == 1
      if(this->mQ == 1 && i_ >= 0)
      {
         auto l1 = this->l();
         ACoeff_t m = n.bottomRows(1) + 1.0;
         // Tau coefficient obtained as
         ACoeff_t f = (2.0*l1 + 2.0*m - 1.0)*(l1 + 2.0*m - 4.0)/(l1 + m - 2.0);
         ACoeff_t nf = (this->normalizeDiag(m, -2, 1)/this->normalizeDiag(m, -2))*f;

         m = n.bottomRows(1) - static_cast<Scalar_t>(k + 1);
         ACoeff_t g;
         switch(k)
         {
            case -1:
               g = this->mI2.d_1(m);
               break;
            case 0:
               g = this->mI2.d0(m);
               break;
            case 1:
               g = this->mI2.d1(m);
               break;
            case 2:
               g = this->mI2.d2(m);
               break;
            default:
               throw std::logic_error("Unknown diagonal for computing correction");
               break;
         }
         ACoeff_t ng = g/this->normalizeDiag(m, k, 1);

         val(i_) -= (nf*ng)(0);
      }
   }

} // Chebyshev
} // Worland
} // SparseSM
} // QuICC
