/**
 * @file SphLapl2Diags.cpp
 * @brief Source of the implementation of the full sphere Bessel SphLapl2Diags
 * sparse operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/Insulating/SphLapl2Diags.hpp"
#include "QuICC/SparseSM/Bessel/Insulating/Utils.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Insulating {

SphLapl2Diags::SphLapl2Diags(const int l) :
    QuICC::SparseSM::Bessel::SphLapl2Diags(BesselKind::INSULATING, l)
{}

SphLapl2Diags::ACoeff_t SphLapl2Diags::d0(const ACoeff_t& n) const
{
   // Compute roots
   std::vector<Scalar_t> roots;
   getRoots(roots, static_cast<int>(this->l()), n.size());

   ACoeff_t val = ACoeff_t::Ones(n.size());
   for (int i = 0; i < roots.size(); i++)
   {
      const auto& k = roots.at(i);
      const auto k2 = k * k;
      val(i) = k2 * k2;
   }

   return val;
}

} // namespace Insulating
} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC
