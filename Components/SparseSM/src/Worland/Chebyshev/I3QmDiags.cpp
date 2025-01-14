/**
 * @file I3QmDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I3QmDiags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I3QmDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   using namespace Internal::Literals;

   I3QmDiags::I3QmDiags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::I3QmDiags(alpha, MHD_MP(-0.5), l, q)
   {
      // q <= 1 is equivalent to no truncation (already zero rows)

      if(q > 1)
      {
         throw std::logic_error("I3Qm: Truncation for q>1 is not implemented");
      }
   }

   I3QmDiags::ACoeff_t I3QmDiags::d_2(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -16.0_mp*(l<1>() + n - 3.0_mp)*(l<1>() + n - 2.0_mp)*(l<1>() + n - 1.0_mp)/((l<1>() + 2.0_mp*n - 5.0_mp)*(l<1>() + 2.0_mp*n - 4.0_mp)*(l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp));

      return this->normalizeDiag(n,-2,-1)*val;
   }

   I3QmDiags::ACoeff_t I3QmDiags::d_1(const ACoeff_t& n) const
   {
      ACoeff_t val;
      val = 8.0_mp*(l<1>() + n - 2.0_mp)*(l<1>() + n - 1.0_mp)*(4.0_mp*l<1>() - 2.0_mp*n - 1.0_mp)/((l<1>() + 2.0_mp*n - 4.0_mp)*(l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp));

      return this->normalizeDiag(n,-1,-1)*val;
   }

   I3QmDiags::ACoeff_t I3QmDiags::d0(const ACoeff_t& n) const
   {
      ACoeff_t val;
      val = -8.0_mp*(l<1>() + n - 1.0_mp)*(2.0_mp*l<2>() - 8.0_mp*l<1>()*n - 2.0_mp*l<1>() - 4.0_mp*n*n + 4.0_mp*n + 3.0_mp)/((l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp));

      return this->normalizeDiag(n,0,-1)*val;
   }

   I3QmDiags::ACoeff_t I3QmDiags::d1(const ACoeff_t& n) const
   {
      ACoeff_t val;
      val = -4.0_mp*(2.0_mp*n + 1.0_mp)*(6.0_mp*l<2>() - 6.0_mp*l<1>() - 4.0_mp*n*n - 4.0_mp*n + 3.0_mp)/((l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp));

      return this->normalizeDiag(n,1,-1)*val;
   }

   I3QmDiags::ACoeff_t I3QmDiags::d2(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -(2.0_mp*n + 1.0_mp)*(2.0_mp*n + 3.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)*(6.0_mp*l<1>() + 2.0_mp*n - 1.0_mp)/((l<1>() + n)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp)*(l<1>() + 2.0_mp*n + 4.0_mp));

      return this->normalizeDiag(n,2,-1)*val;
   }

   I3QmDiags::ACoeff_t I3QmDiags::d3(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -(2.0_mp*n + 1.0_mp)*(2.0_mp*n + 3.0_mp)*(2.0_mp*n + 5.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 3.0_mp)/(2.0_mp*(l<1>() + n)*(l<1>() + n + 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp)*(l<1>() + 2.0_mp*n + 4.0_mp)*(l<1>() + 2.0_mp*n + 5.0_mp));

      return this->normalizeDiag(n,3,-1)*val;
   }

} // Chebyshev
} // Worland
} // SparseSM
} // QuICC
