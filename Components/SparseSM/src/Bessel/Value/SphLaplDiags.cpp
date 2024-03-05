/**
 * @file SphLaplDiags.cpp
 * @brief Source of the implementation of the full sphere Bessel SphLaplDiags
 * sparse operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/Value/SphLaplDiags.hpp"
#include "QuICC/SparseSM/Bessel/Value/Utils.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Value {

SphLaplDiags::SphLaplDiags(const int l) :
    QuICC::SparseSM::Bessel::SphLaplDiags(BesselKind::VALUE, l)
{}

SphLaplDiags::ACoeff_t SphLaplDiags::d0(const ACoeff_t& n) const
{
   // Compute roots
   std::vector<Scalar_t> roots;
   getRoots(roots, static_cast<int>(this->l()), n.size());

   ACoeff_t val = ACoeff_t::Ones(n.size());
   for (int i = 0; i < roots.size(); i++)
   {
      const auto& k = roots.at(i);
      val(i) = -k * k;
   }

   return val;
}

} // namespace Value
} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC
