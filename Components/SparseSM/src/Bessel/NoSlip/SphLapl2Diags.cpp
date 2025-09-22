/**
 * @file SphLapl2Diags.cpp
 * @brief Source of the implementation of the full sphere Bessel SphLapl2Diags
 * sparse operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/NoSlip/SphLapl2Diags.hpp"
#include "QuICC/SparseSM/Bessel/NoSlip/Utils.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace NoSlip {

SphLapl2Diags::SphLapl2Diags(const int l) :
    QuICC::SparseSM::Bessel::SphLapl2Diags(BesselKind::NOSLIP, l)
{}

SphLapl2Diags::ACoeff_t SphLapl2Diags::d0(const ACoeff_t& n) const
{
   // Compute roots
   std::vector<Scalar_t> roots = {0};
   getRoots(roots, static_cast<int>(this->l()), n.size() - 1);

   ACoeff_t val = ACoeff_t::Ones(n.size());
   assert(roots.size() == val.size());
   for (std::size_t i = 0; i < roots.size(); i++)
   {
      const auto& k = roots.at(i);
      const auto k2 = k * k;
      val(i) = k2 * k2;
   }

   return val;
}

} // namespace NoSlip
} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC
