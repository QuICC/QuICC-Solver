/**
 * @file Operators.cpp
 * @brief Source of the Bessel operators
 */

// System include
//
// clang-format off
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"
#include <boost/math/special_functions/bessel.hpp>
// clang-format on

// Project includes
//
#include "QuICC/Polynomial/Bessel/details/Operators.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

namespace details {

Internal::MHDFloat Value_dNu()
{
   using namespace Internal::Literals;
   return 0.5_mp;
}

Internal::MHDFloat Insulating_dNu()
{
   using namespace Internal::Literals;
   return -0.5_mp;
}

Internal::MHDFloat norm(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat dNu)
{
   using namespace Internal::Literals;
   Internal::MHDFloat val;
   unsigned int nu;
   assert(l >= 0);
   if (dNu == Value_dNu())
   {
      nu = l + 1;
   }
   else if (dNu == Insulating_dNu())
   {
      nu = l;
   }
   else
   {
      throw std::logic_error("Unknown Bessel d_nu: " + std::to_string(static_cast<double>(dNu)));
   }
   val = boost::math::sph_bessel(nu, k) / Internal::Math::sqrt(2.0_mp);

   return Internal::Math::abs(val);
}

Internal::MHDFloat SphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   Internal::MHDFloat val;
   val = boost::math::sph_bessel(l, k * r);

   const auto scale = norm(k, l, dNu);
   return val / scale;
}

Internal::MHDFloat rSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   Internal::MHDFloat val;
   val = r * boost::math::sph_bessel(l, k * r);

   const auto scale = norm(k, l, dNu);
   return val / scale;
}

Internal::MHDFloat r_1SphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   using namespace Internal::Literals;
   Internal::MHDFloat val;
   if (l == 0)
   {
      val = boost::math::sph_bessel(l, k * r) / r;
   }
   else
   {
      const Internal::MHDFloat c =
         k / static_cast<Internal::MHDFloat>(2 * l + 1);
      val = c * (boost::math::sph_bessel(l - 1, k * r) +
                   boost::math::sph_bessel(l + 1, k * r));
   }

   const auto scale = norm(k, l, dNu);
   return val / scale;
}

Internal::MHDFloat dSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   auto dl = static_cast<Internal::MHDFloat>(l);
   Internal::MHDFloat val;
   if (l == 0)
   {
      val = (dl / r) * boost::math::sph_bessel(l, k * r) -
            k * boost::math::sph_bessel(l + 1, k * r);
   }
   else
   {
      const Internal::MHDFloat c =
         k / static_cast<Internal::MHDFloat>(2 * l + 1);
      val = (c * dl) * boost::math::sph_bessel(l - 1, k * r) +
            (c * dl - k) * boost::math::sph_bessel(l + 1, k * r);
   }

   const auto scale = norm(k, l, dNu);
   return val / scale;
}

Internal::MHDFloat drSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   auto dl1 = static_cast<Internal::MHDFloat>(l + 1);
   Internal::MHDFloat val;
   val = dl1 * boost::math::sph_bessel(l, k * r) -
         k * r * boost::math::sph_bessel(l + 1, k * r);

   const auto scale = norm(k, l, dNu);
   return val / scale;
}

Internal::MHDFloat r_1drSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   auto dl1 = static_cast<Internal::MHDFloat>(l + 1);
   Internal::MHDFloat val;
   if (l == 0)
   {
      val = dl1 * boost::math::sph_bessel(l, k * r) / r -
            k * boost::math::sph_bessel(l + 1, k * r);
   }
   else
   {
      const Internal::MHDFloat c =
         k / static_cast<Internal::MHDFloat>(2 * l + 1);
      val = dl1 * c * boost::math::sph_bessel(l - 1, k * r) +
            (dl1 * c - k) * boost::math::sph_bessel(l + 1, k * r);
   }

   const auto scale = norm(k, l, dNu);
   return val / scale;
}

Internal::MHDFloat slaplSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   Internal::MHDFloat val;
   val = -k * k * boost::math::sph_bessel(l, k * r);

   const auto scale = norm(k, l, dNu);
   return val / scale;
}

Internal::MHDFloat raiseSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   using namespace Internal::Literals;
   Internal::MHDFloat val;

   val = k * boost::math::sph_bessel(l + 1, k * r);

   const auto scale = norm(k, l, dNu);
   return val / scale;
}

Internal::MHDFloat lowerSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   using namespace Internal::Literals;
   Internal::MHDFloat val;

   if (l > 0)
   {
      val = k * boost::math::sph_bessel(l - 1, k * r);
   }
   else
   {
      throw std::logic_error("Negative Bessel parameter is not defined");
   }

   const auto scale = norm(k, l, dNu);
   return val / scale;
}

void getRoots(std::vector<Internal::MHDFloat>& roots, const int l,
   const int nRoots, const Internal::MHDFloat dNu)
{
   Internal::MHDFloat nu = static_cast<Internal::MHDFloat>(l) + dNu;
   boost::math::cyl_bessel_j_zero(nu, 1, nRoots, std::back_inserter(roots));
}

} // namespace details
} // namespace Bessel
} // namespace Polynomial
} // namespace QuICC
