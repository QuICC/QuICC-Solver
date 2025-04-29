/**
 * @file SphLaplDiags.cpp
 * @brief Source of the implementation of the full sphere Bessel SphLaplDiags
 * sparse operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/Insulating/SphLaplDiags.hpp"
#include "QuICC/SparseSM/Bessel/Insulating/Utils.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Insulating {

SphLaplDiags::SphLaplDiags(const int l) :
    QuICC::SparseSM::Bessel::SphLaplDiags(BesselKind::INSULATING, l)
{}

SphLaplDiags::ACoeff_t SphLaplDiags::d0(const ACoeff_t& n) const
{
   // Compute roots
   std::vector<Scalar_t> roots;
   getRoots(roots, static_cast<int>(this->l()), n.size());

   ACoeff_t val = ACoeff_t::Ones(n.size());
   for (std::size_t i = 0; i < roots.size(); i++)
   {
      const auto& k = roots.at(i);
      val(i) = -k * k;
   }

   return val;
}

} // namespace Insulating
} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC
