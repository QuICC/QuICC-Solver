/**
 * @file JacobiBase.hpp
 * @brief Implementation of the Jacobi polynomial base
 */

#ifndef QUICC_POLYNOMIAL_JACOBI_JACOBIASYMPTOTIC_HPP
#define QUICC_POLYNOMIAL_JACOBI_JACOBIASYMPTOTIC_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Polynomial/Jacobi/JacobiBase.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Polynomial {

namespace Jacobi {

namespace JacobiAsy {

/**
 * @brief boundary asymptotic evaluation factor for weights eq 3.24 [1]
 *
 * Refs
 *
 * [1] Fast and accurate computation of Gauss–Legendre and Gauss–Jacobi quadrature nodes and weights
 *     Nicholas Hale and Alex Townsend 2013
 * [2] The bounds for the error term of an asymptotic approximation of Jacobi polynomials
 *     Barella e Gatteschi 1998
 *
 */
std::array<MHDFloat, 2> evalPdPatBnd(const std::uint32_t n, const MHDFloat a, const MHDFloat b, const MHDFloat t);

namespace details {

/**
 * @brief compute g, g', g''
 *
 */
inline std::array<MHDFloat, 3> g(const MHDFloat a, const MHDFloat b, const MHDFloat t)
{
   constexpr MHDFloat tLimit = 0.4;
   // g(theta)
   auto gc1 = 0.25 - a*a;
   auto gc2 = 0.25 - b*b;
   auto t_2 = static_cast<MHDLong>(0.5*t);
   auto tan_t_2 = std::tan(t_2);
   auto cot_t_2 = 1.0/tan_t_2;
   auto tP2 = static_cast<MHDLong>(t)*static_cast<MHDLong>(t);
   auto tP4 = tP2*tP2;
   auto tP8 = tP4*tP4;
   auto tG = cot_t_2 - 2.0/t;
   // using Taylor series for small t
   if (t < tLimit)
   {
      tG = - t/6. - tP2*t/360. -tP4*t/15120. -tP4*tP2*t/604800.
         - tP8*t/23950080. - tP8*tP2*t*691./653837184000.
         - tP8*tP4*t/37362124800.;
   }
   auto g = gc1*tG - gc2*tan_t_2;
   // g'
   auto csc_t_2 = 1.0/std::sin(t_2);
   auto sec_t_2 = 1.0/std::cos(t_2);
   auto csc_t_2P2 = csc_t_2*csc_t_2;
   auto sec_t_2P2 = sec_t_2*sec_t_2;
   // using Taylor series for small t
   auto tGp = 2.0/(tP2) - 0.5*csc_t_2P2;
   if (t < tLimit)
   {
      tGp = -(1./6. + tP2/120. + tP4/3024. + tP4*tP2/86400.
         + tP8/2661120. + 691.*tP8*tP2/59439744000.
         + tP8*tP4/2874009600 + 3617.*tP8*tP4*tP2/355687428096000.
         + 43867*tP8*tP8/150267476975616000);
   }
   auto gp = gc1*tGp - 0.5*gc2*sec_t_2P2;
   // g''
   // using Taylor series for small t
   auto tGpp = -4.0/(tP2*t) + 0.5*cot_t_2*csc_t_2P2;
   if (t < tLimit)
   {
      tGpp = -t/60. - tP2*t/756. -tP4*t/14400. -tP4*tP2*t/332640.
         -691.*tP8*t/5943974400. -tP8*tP2*t/239500800.
         -41056567.*tP8*tP4*t/285665824604160000.
         -20728717.*tP8*tP4*tP2*t/4627786358587392000.;
   }
   auto gpp = gc1*tGpp
      - 0.5*gc2*sec_t_2P2*tan_t_2;

   return std::array<MHDFloat, 3>{static_cast<MHDFloat>(g), static_cast<MHDFloat>(gp), static_cast<MHDFloat>(gpp)};
}

/**
 * @brief return Stirling's coefficients
 *
 * numerator https://oeis.org/A001163
 * denominator https://oeis.org/A001164
 *
 */
constexpr inline std::array<MHDFloat, 10> stirling()
{
   return std::array<MHDFloat, 10>
   {
      1., 1./12., 1./288.,
      -139./51840., -571./2488320., 163879./209018880.,
      5246819./75246796800., -534703531./902961561600., -4483131259./86684309913600.,
      432261921612371./514904800886784000.
   };
}

inline MHDFloat expArgBnd(const std::uint32_t n, const MHDFloat a)
{
   // Exponent argument
   auto epsilon = std::numeric_limits<decltype(a)>::epsilon();
   static constexpr int ulp = 2;
   constexpr std::uint32_t M = 50;
   auto expNum = a*a;
   auto expDen = static_cast<decltype(a)>(n);
   auto expArg = expNum/(2.0*expDen);
   auto sgn = 1.0;
   for (std::uint32_t j=2; j<M; ++j)
   {
      expNum *= a;
      expDen *= n;
      sgn = -sgn;
      auto delta = sgn*expNum/(expDen*j*(j+1));
      expArg += delta;

      auto tol = epsilon * std::abs(expArg) * ulp;
      if (std::abs(delta) <= tol)
      {
         break;
      }
   }
   return expArg;
}

inline MHDFloat expArgInt(const std::uint32_t n, const MHDFloat a, const MHDFloat b)
{
   // Exponent argument
   auto epsilon = std::numeric_limits<decltype(a)>::epsilon();
   static constexpr int ulp = 2;
   constexpr std::uint32_t M = 50;
   auto one_n = 1.0/static_cast<decltype(a)>(n);
   auto coeffN = -1.0;
   auto coeffQ = 1.0;
   auto deltaA = a;
   auto deltaB = b;
   auto deltaAB = a+b;
   auto expArg = 0.0;
   for (std::uint32_t j=1; j<M; ++j)
   {
      coeffQ = j*(j+1);
      coeffN *= -one_n;
      deltaA *= a;
      deltaB *= b;
      deltaAB *= (a+b)/2.0;
      auto delta = coeffN/coeffQ*(deltaA + deltaB - deltaAB);
      expArg += delta;

      auto tol = epsilon * std::abs(expArg) * ulp;
      // always end on even iteration to avoid cancellation if a = -b
      if (std::abs(delta) <= tol && (j % 2 != 0))
      {
         break;
      }
   }
   return expArg;
}


/**
 * @brief compute f(theta) eq 3.23 [1]
 *
 */
inline MHDFloat f(const std::uint32_t n, const MHDFloat a, const MHDFloat b, const MHDFloat t, const std::uint32_t m)
{
   auto rho = n + 0.5*(a+b+1);
   auto phiM = 0.5*(2*rho+m)*t;

   auto t_2 = 0.5*t;
   auto sint_2 = std::sin(t_2);
   auto cost_2 = std::cos(t_2);

   auto f = 0.0;
   for (std::uint32_t l = 0; l <= m; ++l)
   {

      // pow/fact/Poch l
      auto sint_2l = 1.0;
      auto lf = 1.0;
      auto aPlus_poch = 1.0;
      auto aMinus_poch = 1.0;
      for (std::uint32_t i = 1; i <= l; ++i)
      {
         lf *= i;
         sint_2l *= sint_2;
         aPlus_poch *= 0.5+a+i-1;
         aMinus_poch *= 0.5-a+i-1;
      }

      // pow/fact ml
      auto cost_2ml = 1.0;
      auto mlf = 1.0;
      auto bPlus_poch = 1.0;
      auto bMinus_poch = 1.0;
      for (std::uint32_t i = 1; i <= m-l; ++i)
      {
         mlf *= i;
         cost_2ml *= cost_2;
         bPlus_poch *= 0.5+b+i-1;
         bMinus_poch *= 0.5-b+i-1;
      }

      auto c = aPlus_poch*aMinus_poch*bPlus_poch*bMinus_poch;
      auto phi = phiM - 0.5*(a+l+0.5)*Math::PI;
      auto cosphi = std::cos(phi);
      f+= c*cosphi/(lf*mlf*sint_2l*cost_2ml);

   }
   return f;
}

} // namespace details

/**
 * @brief interior asymptotic evaluation factor for weights eq 3.22 [1]
 *
 */
std::array<MHDFloat, 2> evalPdPatInt(const std::uint32_t n, const MHDFloat a, const MHDFloat b, const MHDFloat t);

} // namespace JacobiAsy
} // namespace Jacobi
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_JACOBI_JACOBIASYMPTOTIC_HPP
