/**
 * @file R1D1DivR1.cpp
 * @brief Source of the implementation of boundary value for toroidal stress-free
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/Boundary/R1D1DivR1.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Typedefs.hpp"
#include "Polynomial/SphericalBessel/Jnl.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Boundary {

   R1D1DivR1::R1D1DivR1(const BesselKind type, const int l)
      : ICondition(type, l)
   {
   }

   R1D1DivR1::ACoeff_t R1D1DivR1::compute(const int maxN)
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
         const auto& k = roots.at(i);
         val(i) = Polynomial::SphericalBessel::SphJnl(k, il, 1_mp, dNu) - 
            Polynomial::SphericalBessel::dSphJnl(k, il, 1_mp, dNu);
      }

      return val;
   }

} // Boundary
} // Bessel
} // SparseSM
} // QuICC
