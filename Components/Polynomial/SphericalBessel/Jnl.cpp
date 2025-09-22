/**
 * @file Jnl.cpp
 * @brief Source of the spherical Bessel operators
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
#include "Polynomial/SphericalBessel/Jnl.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace Polynomial {

namespace SphericalBessel {

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

Internal::MHDFloat NoSlip_dNu()
{
   using namespace Internal::Literals;
   return 1.5_mp;
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
   else if (dNu == NoSlip_dNu())
   {
      nu = l + 2;
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

   if(k == 0)
   {
      val = Fundamental(l, r);
   }
   else
   {
      val = boost::math::sph_bessel(l, k * r)/norm(k, l, dNu);
   }

   return val;
}

Internal::MHDFloat rSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   Internal::MHDFloat val;

   if(k == 0)
   {
      val = rFundamental(l, r);
   }
   else
   {
      val = r * boost::math::sph_bessel(l, k * r)/norm(k, l, dNu);
   }
   return val;
}

Internal::MHDFloat r_1SphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   Internal::MHDFloat val;
   if(k == 0)
   {
      val = r_1Fundamental(l, r);
   }
   else
   {
      if (l == 0)
      {
         val = boost::math::sph_bessel(l, k * r) / r;
      }
      else
      {
         using namespace Internal::Literals;
         const Internal::MHDFloat c =
            k / static_cast<Internal::MHDFloat>(2 * l + 1);
         val = c * (boost::math::sph_bessel(l - 1, k * r) +
                      boost::math::sph_bessel(l + 1, k * r));
      }

      val /= norm(k, l, dNu);
   }

   return val;
}

Internal::MHDFloat dSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   Internal::MHDFloat val;
   if(k == 0)
   {
      val = dFundamental(l, r);
   }
   else
   {
      auto dl = static_cast<Internal::MHDFloat>(l);
      if (l == 0)
      {
         val = - k * boost::math::sph_bessel(1, k * r);
      }
      else
      {
         const Internal::MHDFloat c =
            k / static_cast<Internal::MHDFloat>(2 * l + 1);
         val = (c * dl) * boost::math::sph_bessel(l - 1, k * r) +
               (c * dl - k) * boost::math::sph_bessel(l + 1, k * r);
      }

      val /= norm(k, l, dNu);
   }

   return val;
}

Internal::MHDFloat d2SphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   Internal::MHDFloat val;
   if(k == 0)
   {
      val = d2Fundamental(l, r);
   }
   else
   {
      using namespace Internal::Literals;
      if (l == 0)
      {
         auto c = - k * k / 3_mp;
         val = c * (boost::math::sph_bessel(0, k * r) -
               2_mp*boost::math::sph_bessel(2, k * r));
      }
      else if (l == 1)
      {
         auto c = k * k / 5_mp;
         val = c * (- 3_mp * boost::math::sph_bessel(1, k * r) +
               2_mp*boost::math::sph_bessel(3, k * r));
      }
      else
      {
         auto dl = static_cast<Internal::MHDFloat>(l);
         const Internal::MHDFloat c =
            k*k / ((2_mp*dl - 1_mp)*(2_mp*dl + 1_mp)*(2_mp*dl + 3_mp));
         val = c*(
               (dl - 1_mp)*dl*(2_mp*dl + 3_mp)*boost::math::sph_bessel(l - 2, k * r) -
               (2_mp*dl*dl + 2_mp*dl - 1_mp)*(2_mp*dl + 1_mp)*boost::math::sph_bessel(l, k * r) +
               (dl+1_mp)*(dl+2_mp)*(2_mp*dl-1_mp)*boost::math::sph_bessel(l + 2, k * r)
               );
      }

      val /= norm(k, l, dNu);
   }

   return val;
}

Internal::MHDFloat drSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   Internal::MHDFloat val;

   if(k == 0)
   {
      val = drFundamental(l, r);
   }
   else
   {
      auto dl1 = static_cast<Internal::MHDFloat>(l + 1);
      val = dl1 * boost::math::sph_bessel(l, k * r) -
            k * r * boost::math::sph_bessel(l + 1, k * r);

      val /= norm(k, l, dNu);
   }

   return val;
}

Internal::MHDFloat r_1drSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   Internal::MHDFloat val;

   if(k == 0)
   {
      val = r_1drFundamental(l, r);
   }
   else
   {
      auto dl = static_cast<Internal::MHDFloat>(l);
      auto dl1 = static_cast<Internal::MHDFloat>(l + 1);
      if (l == 0)
      {
         val = dl1 * boost::math::sph_bessel(l, k * r) / r -
               k * boost::math::sph_bessel(l + 1, k * r);
      }
      else
      {
         const Internal::MHDFloat c =
            k / static_cast<Internal::MHDFloat>(2 * l + 1);

         val = c*(dl1 * boost::math::sph_bessel(l - 1, k * r) -
               dl * boost::math::sph_bessel(l + 1, k * r));
      }

      val /= norm(k, l, dNu);
   }

   return val;
}

Internal::MHDFloat slaplSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu)
{
   Internal::MHDFloat val;

   if(k == 0)
   {
      val = slaplFundamental(l, r);
   }
   else
   {
      val = -k * k * boost::math::sph_bessel(l, k * r);

      val /= norm(k, l, dNu);
   }

   return val;
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

Internal::MHDFloat normFundamental(const int l)
{
   using namespace Internal::Literals;
   Internal::MHDFloat val;
   assert(l >= 0);
   val = 1_mp/Internal::Math::sqrt(3_mp + 2_mp*l);

   return val;
}

Internal::MHDFloat Fundamental(const int l, const Internal::MHDFloat r)
{
   Internal::MHDFloat val;

   val = Internal::Math::pow(r,l);
   val /= normFundamental(l);

   return val;
}

Internal::MHDFloat rFundamental(const int l, const Internal::MHDFloat r)
{
   Internal::MHDFloat val;

   val = Internal::Math::pow(r,l+1);
   val /= normFundamental(l);

   return val;
}

Internal::MHDFloat r_1Fundamental(const int l, const Internal::MHDFloat r)
{
   Internal::MHDFloat val;

   val = Internal::Math::pow(r,l-1);
   val /= normFundamental(l);

   return val;
}

Internal::MHDFloat dFundamental(const int l, const Internal::MHDFloat r)
{
   Internal::MHDFloat val;

   if(l < 1)
   {
      val = 0;
   }
   else
   {
      val = static_cast<Internal::MHDFloat>(l)*Internal::Math::pow(r,l-1);
      val /= normFundamental(l);
   }

   return val;
}

Internal::MHDFloat d2Fundamental(const int l, const Internal::MHDFloat r)
{
   Internal::MHDFloat val;

   if(l < 2)
   {
      val = 0;
   }
   else
   {
      val = static_cast<Internal::MHDFloat>(l*(l-1))*Internal::Math::pow(r,l-2);
      val /= normFundamental(l);
   }

   return val;
}

Internal::MHDFloat drFundamental(const int l, const Internal::MHDFloat r)
{
   Internal::MHDFloat val;

   val = static_cast<Internal::MHDFloat>(l+1)*Internal::Math::pow(r,l);
   val /= normFundamental(l);

   return val;
}

Internal::MHDFloat r_1drFundamental(const int l, const Internal::MHDFloat r)
{
   Internal::MHDFloat val;

   val = static_cast<Internal::MHDFloat>(l+1)*Internal::Math::pow(r,l-1);
   val /= normFundamental(l);

   return val;
}

Internal::MHDFloat slaplFundamental(const int l, const Internal::MHDFloat r)
{
   Internal::MHDFloat val = 0;

   return val;
}

} // namespace SphericalBessel
} // namespace Polynomial
} // namespace QuICC
