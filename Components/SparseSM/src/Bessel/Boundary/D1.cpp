/**
 * @file D1.cpp
 * @brief Source of the implementation of boundary value of first derivative
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/Boundary/D1.hpp"
#include "Polynomial/SphericalBessel/Jnl.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Boundary {

D1::D1(const BesselKind type, const int l) : ICondition(type, l) {}

D1::ACoeff_t D1::compute(const int maxN)
{
   using namespace Internal::Literals;

   int nRoot = maxN + 1;
   ACoeff_t val = ACoeff_t::Ones(nRoot);
   Internal::MHDFloat dNu;
   std::vector<Internal::MHDFloat> roots;
   if (this->type() == BesselKind::VALUE)
   {
      dNu = Polynomial::SphericalBessel::Value_dNu();
   }
   else if (this->type() == BesselKind::INSULATING)
   {
      dNu = Polynomial::SphericalBessel::Insulating_dNu();
   }
   else if (this->type() == BesselKind::NOSLIP)
   {
      dNu = Polynomial::SphericalBessel::NoSlip_dNu();
      roots.push_back(0);
      nRoot--;
   }
   else
   {
      throw std::logic_error("Unknown Bessel Kind for D1 boundary");
   }

   int il = static_cast<int>(this->l());
   Polynomial::SphericalBessel::getRoots(roots, il, nRoot, dNu);
   for (int i = 0; i < val.size(); i++)
   {
      const auto& k = roots.at(i);
      val(i) = Polynomial::SphericalBessel::dSphJnl(k, il, 1_mp, dNu);
   }

   return val;
}

} // namespace Boundary
} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC
