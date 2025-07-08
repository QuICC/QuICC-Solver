/**
 * @file InsulatingSphere.cpp
 * @brief Source of the implementation of boundary value for an insulating outer sphere
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/Boundary/InsulatingSphere.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Typedefs.hpp"
#include "Polynomial/SphericalBessel/Jnl.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Boundary {

   InsulatingSphere::InsulatingSphere(const BesselKind type, const int l)
      : ICondition(type, l)
   {
   }

   InsulatingSphere::ACoeff_t InsulatingSphere::compute(const int maxN)
   {
      using namespace Internal::Literals;

      ACoeff_t val = ACoeff_t::Ones(maxN+1);
      Internal::MHDFloat dNu;
      if(this->type() == BesselKind::VALUE)
      {
         dNu = Polynomial::SphericalBessel::Value_dNu();
      }
      else if(this->type() == BesselKind::INSULATING)
      {
         dNu = Polynomial::SphericalBessel::Insulating_dNu();
      }
      else
      {
         throw std::logic_error("Unknown Bessel Kind");
      }

      std::vector<Internal::MHDFloat> roots;
      int il = static_cast<int>(this->l());
      Polynomial::SphericalBessel::getRoots(roots, il, maxN+1, dNu);
      for(int i = 0; i < val.size(); i++)
      {
         auto l1 = this->l() + 1_mp;
         const auto& k = roots.at(i);
         val(i) = Polynomial::SphericalBessel::dSphJnl(k, il, 1_mp, dNu) + 
         + l1*Polynomial::SphericalBessel::SphJnl(k, il, 1_mp, dNu);
      }

      return val;
   }

} // Boundary
} // Bessel
} // SparseSM
} // QuICC
