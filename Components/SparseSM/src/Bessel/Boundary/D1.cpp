/**
 * @file D1.cpp
 * @brief Source of the implementation of boundary value of first derivative
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/Boundary/D1.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Typedefs.hpp"
#include "Polynomial/SphericalBessel/Jnl.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Boundary {

   D1::D1(const BesselKind type, const int l)
      : ICondition(type, l)
   {
   }

   D1::ACoeff_t D1::compute(const int maxN)
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
      Polynomial::SphericalBessel::getRoots(roots, this->l(), maxN+1, dNu);
      for(int i = 0; i < val.size(); i++)
      {
         const auto& k = roots.at(i);
         val(i) = Polynomial::SphericalBessel::dSphJnl(k, this->l(), 1_mp, dNu);
      }

      return val;
   }

} // Boundary
} // Bessel
} // SparseSM
} // QuICC
