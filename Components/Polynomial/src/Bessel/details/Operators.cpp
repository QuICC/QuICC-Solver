/**
 * @file Operators.cpp
 * @brief Source of the Bessel operators
 */

// System include
//
#include <boost/math/special_functions/bessel.hpp>

// Project includes
//
#include "QuICC/Polynomial/Bessel/details/Operators.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

namespace details {

   Internal::MHDFloat SphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r)
   {
      Internal::MHDFloat val;
      val = boost::math::sph_bessel(l, k*r);

      return val;
   }

   Internal::MHDFloat rSphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r)
   {
      Internal::MHDFloat val;
      val = r*boost::math::sph_bessel(l, k*r);

      return val;
   }

   Internal::MHDFloat r_1SphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r)
   {
      Internal::MHDFloat val;
      val = boost::math::sph_bessel(l, k*r)/r;

      return val;
   }

   Internal::MHDFloat dSphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r)
   {
      auto dl = static_cast<Internal::MHDFloat>(l);
      Internal::MHDFloat val;
      val = (dl/r)*boost::math::sph_bessel(l, k*r) - k*boost::math::sph_bessel(l + 1, k*r);

      return val;
   }

   Internal::MHDFloat drSphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r)
   {
      auto dl1 = static_cast<Internal::MHDFloat>(l+1);
      Internal::MHDFloat val;
      val = dl1 * boost::math::sph_bessel(l, k * r) - k * r * boost::math::sph_bessel(l + 1, k * r);

      return val;
   }

   Internal::MHDFloat r_1drSphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r)
   {
      auto dl1 = static_cast<Internal::MHDFloat>(l+1);
      Internal::MHDFloat val;
      val = dl1 * boost::math::sph_bessel(l, k * r)/r - k * boost::math::sph_bessel(l + 1, k * r);

      return val;
   }

   Internal::MHDFloat slaplSphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r)
   {
      Internal::MHDFloat val;
      val = -k*k*boost::math::sph_bessel(l, k*r);

      return val;
   }

   void getTorRoots(std::vector<Internal::MHDFloat>& roots, const int l, const int nRoots)
   {
      using namespace Internal::Literals;
      Internal::MHDFloat nu = static_cast<Internal::MHDFloat>(l) + 0.5_mp;
      boost::math::cyl_bessel_j_zero(nu, 1, nRoots, std::back_inserter(roots));
   }

   void getPolRoots(std::vector<Internal::MHDFloat>& roots, const int l, const int nRoots)
   {
      using namespace Internal::Literals;
      Internal::MHDFloat nu = static_cast<Internal::MHDFloat>(l) - 0.5_mp;
      boost::math::cyl_bessel_j_zero(nu, 1, nRoots, std::back_inserter(roots));
   }

}
}
}
}
