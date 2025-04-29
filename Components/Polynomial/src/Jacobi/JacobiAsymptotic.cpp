/**
 * @file JacobiAsymptotic.cpp
 * @brief Implementation of the Jacobi polynomial base
 */

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
#include "QuICC/Polynomial/Jacobi/JacobiAsymptotic.hpp"
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"
#include "Bessel/Jn.hpp"


namespace QuICC {

namespace Polynomial {

namespace Jacobi {

namespace JacobiAsy {


std::array<MHDFloat, 2> evalPdPatBnd(const std::uint32_t n, const MHDFloat a, const MHDFloat b, const MHDFloat t)
{
   // preconditions
   assert(a >= -1.0);
   assert(b >= -1.0);
   assert(t >= 0);

   // g(theta), g', g''
   std::array<MHDFloat, 3> gs = details::g(a, b, t);
   auto g = gs[0];
   auto gp = gs[1];

   // sums A_m and B_m [2]
   auto A0 = 1.0;
   auto tB0 = 0.25*g;
   auto gc1 = 0.25 - a*a;
   auto gc2 = 0.25 - b*b;
   auto A1at0 = a*(gc1+3.*gc2)/24.;
   auto onePlusTwoA = 1.0+2.0*a;
   auto A1 = 0.125*(gp - onePlusTwoA*g/t -0.25*g*g) - A1at0;

   // higher order terms require numerical diff and int of A_m(theta), B_m(theta)

   // get Chebyshev GL points [0, theta]
   constexpr std::uint32_t nPts = 8;
   Array grid(nPts), weights(nPts);
   Quadrature::ChebyshevRule<Quadrature::gaussLobatto_t> quad;
   quad.computeQuadrature(grid, weights, nPts, t, 0);
   // avoid division by zero
   grid(0) = std::numeric_limits<Internal::MHDFloat>::epsilon()*t;
   Matrix D(nPts, nPts), I(nPts, nPts), T(nPts, nPts), TI(nPts, nPts);

   // collocation differentiation matrix
   quad.computeDiffMat(D);
   D *= 2.0/t; // scale to domain
   // integration matrix
   quad.computeIntMat(I);
   // forward projection matrix (modal -> phys)
   quad.computeTMat<1>(T);
   // backward projecton matrix (phys -> modal)
   quad.computeTMat<0>(TI);
   // collocation integration matrix
   auto Q = 0.5*t*T*I*TI;


   // eval A1'/t and f*A1 over [0,theta]
   Array A1p(nPts), A1p_t(nPts), f(nPts), fmulA1(nPts);
   for (std::size_t i = 0; i < nPts; ++i)
   {
      auto ti = grid(i);
      std::array<MHDFloat, 3> gs = details::g(a, b, ti);
      // A1'
      A1p(i) = 0.125*gs[2] + onePlusTwoA*(gs[0]/(8.0*ti*ti) - gs[1]/(8.0*ti))
      - 0.0625*gs[0]*gs[1];
      // A1'/t
      A1p_t(i) = (0.125*gs[2] + onePlusTwoA*(gs[0]/(8.0*ti*ti) - gs[1]/(8.0*ti))
         - 0.0625*gs[0]*gs[1] )/ti;
      // f = 1/2 * gp
      f(i) = 0.5*gs[1];
      // f*A1 = 1/2 * gp * A1
      fmulA1(i) = 0.5*gs[1]
         * (0.125*(gs[1] - onePlusTwoA*gs[0]/ti -0.25*gs[0]*gs[0]) - A1at0);
   }

   // limit t->0
   A1p_t(0) = -gc1/720. - gc1*gc1/576. - gc1*gc2/96. - gc2*gc2/64. - gc2/48.
      + a*(gc1/720. + gc2/48.);
   auto IA1p_t = Q*A1p_t;
   auto IfmulA1 = Q*fmulA1;

   // B1
   auto halfPlusA = 0.5+a;
   auto tB1 = (-0.5*A1p - halfPlusA*IA1p_t + 0.5*IfmulA1).eval();
   tB1(0) = 0;
   auto B1 = (tB1.array()/grid.array()).matrix().eval();
   // limit t->0
   B1(0) =  gc1/720. + gc1*gc1/576. + gc1*gc2/96. + gc2*gc2/64. + gc2/48.
      + a*(gc1*gc1/576. + gc2*gc2/64. + gc1*gc2/96.)
      - a*a*(gc1/720. + gc2/48.);

   // A2
   auto IfmultB1 = Q*(f.array()*tB1.array()).matrix();
   auto DtB1 = D*tB1;
   auto A2 = (0.5*DtB1 - halfPlusA*B1 -0.5*IfmultB1).eval();
   // constrain
   A2 = A2.array() - A2(0);

   // A2p
   auto A2p = (D*A2).eval();
   A2p = A2p.array() - A2p(0);
   auto A2p_t = (A2p.array()/grid.array()).matrix().eval();
   auto nPtsM1 = nPts-1;
   // extrapolate t = 0
   {
      auto nom = 0.0;
      auto den = 0.0;
      auto sgn = 1.0;
      for (std::size_t i = 1; i < nPts; ++i)
      {
         auto w = sgn*(Math::PI*0.5-grid(i));
         if (i == nPtsM1)
         {
            w *= 0.5;
         }
         nom += w*A2p_t(i);
         den += w;
         sgn = -sgn;
      }
      A2p_t(0) = nom/den;
   }

   // B2
   auto tB2 = - 0.5*A2p - halfPlusA*Q*A2p_t + 0.5*Q*(f.array()*A2.array()).matrix();
   auto B2 = (tB2.array()/grid.array()).matrix().eval();
   // extrapolate at t = 0
   {
      auto nom = 0.0;
      auto den = 0.0;
      auto sgn = 1.0;
      for (std::size_t i = 1; i < nPts; ++i)
      {
         auto w = sgn*(Math::PI*0.5-grid(i));
         if (i == nPtsM1)
         {
            w *= 0.5;
         }
         nom += w*B2(i);
         den += w;
         sgn = -sgn;
      }
      B2(0) = nom/den;
   }

   // A3
   auto IfmultB2 = Q*(f.array()*tB2.array()).matrix();
   auto DtB2 = D*tB2;
   auto A3 = (0.5*DtB2 - halfPlusA*B2 - 0.5*IfmultB2).eval();
   // constrain
   A3 = A3.array() - A3(0);

   auto rho = n + 0.5*(a + b + 1.0);
   auto rhoP2 = rho*rho;
   auto sumA = A0 + A1/(rhoP2) + A2(nPtsM1)/(rhoP2*rhoP2)
      + A3(nPtsM1)/(rhoP2*rhoP2*rhoP2);
   auto sumtB = tB0/rho + tB1(nPtsM1)/(rhoP2*rho)
      + tB2(nPtsM1)/(rhoP2*rhoP2*rho);

   auto t_2 = 0.5*t;
   auto lhs = std::pow(std::sin(t_2), a+0.5)*std::pow(std::cos(t_2), b+0.5);

   // Stirling's ratio
   auto sna = 1.0;
   auto sn = 1.0;
   auto powNa = 1.0;
   auto powN = 1.0;
   auto strlng = details::stirling();
   for (std::size_t i = 1; i < 10; ++i)
   {
      powNa *= n+a;
      powN *= n;
      sna += strlng[i]/powNa;
      sn += strlng[i]/powN;
   }

   // exponent argument
   auto expArgBnd = details::expArgBnd(n, a);
   // prefactor on rhs
   auto rhs = std::sqrt((n+a)*t/(n*2.0))*std::exp(expArgBnd)*std::pow(n/rho, a)*sna/sn;
   // Bessel
   using namespace QuICC::Polynomial::Bessel;
   auto Ja = besselJ(a, rho*t);
   auto Jb = besselJ(a+1, rho*t);
   // equation (3.24) [1]
   auto P = rhs*(Ja*sumA + Jb*sumtB) / lhs;

   // P_{n-1}
   // sums A_m and B_m [2]
   auto rho_m1 = rho-1.0;
   auto rho_m1P2 = rho_m1*rho_m1;
   auto sumA_m1 = A0 + A1/rho_m1P2 + A2(nPtsM1)/(rho_m1P2*rho_m1P2)
      + A3(nPtsM1)/(rho_m1P2*rho_m1P2*rho_m1P2);
   auto sumtB_m1 = tB0/rho_m1 + tB1(nPtsM1)/(rho_m1P2*rho_m1)
      + tB2(nPtsM1)/(rho_m1P2*rho_m1P2*rho_m1);
   // missing higher order terms...
   // prefactor on rhs
   // instead of recomputing the exp arg, we flip the term in the sqrt
   auto rhs_m1 = std::sqrt(n*t/((n+a)*2.0))*std::exp(expArgBnd)*std::pow(n/rho_m1, a)*sna/sn;
   // Bessel
   auto Ja_m1 = besselJ(a, rho_m1*t);
   auto Jb_m1 = besselJ(a+1, rho_m1*t);
   // equation (3.24) [1]
   auto P_m1 = rhs_m1*(Ja_m1*sumA_m1 + Jb_m1*sumtB_m1) / lhs;

   // first derivative
   auto cs = JacobiBase::dPnabPnab<Quadrature::natural_t>(n, a, b);
   auto dP = ((cs(1)+cs(2)*t)*P + cs(3)*P_m1)/(-cs(0)*std::sin(t));

   return std::array<MHDFloat, 2>{P, dP};
}


std::array<MHDFloat, 2> evalPdPatInt(const std::uint32_t n, const MHDFloat a, const MHDFloat b, const MHDFloat t)
{
   // preconditions
   assert(a >= -1.0);
   assert(b >= -1.0);
   assert(t > 0);

   auto rho = n + 0.5*(a + b + 1.0);
   auto rho_m1 = rho - 1.0;

   // Sum f
   auto epsilon = std::numeric_limits<decltype(a)>::epsilon();
   static constexpr int ulp = 2;
   constexpr std::uint32_t M = 20;
   auto sumF = 0.0;
   auto sumF_m1 = 0.0;
   auto one_two_pow = 1.0;
   auto twoRhoPlusone_poch = 1.0;
   auto twoRho_m1Plusone_poch = 1.0;
   for (std::size_t i = 0; i < M; ++i)
   {
      auto delta = details::f(n, a, b, t, i) * one_two_pow / twoRhoPlusone_poch;
      sumF += delta;

      auto delta_m1 = details::f(n-1, a, b, t, i) * one_two_pow / twoRho_m1Plusone_poch;
      sumF_m1 += delta_m1;


      // Stop if converged
      auto tol = epsilon * std::abs(sumF) * ulp;
      if (i > 10 && std::abs(delta) <= tol)
      {
         break;
      }

      // update exp/Poch
      one_two_pow *= 0.5;
      twoRhoPlusone_poch *= 2.0*rho+1.0 + i;
      twoRho_m1Plusone_poch *= 2.0*rho_m1+1.0 + i;
   }

   // Stirling's ratio
   auto sna = 1.0;
   auto snb = 1.0;
   auto s2nab = 1.0;
   auto powNa = 1.0;
   auto powNb = 1.0;
   auto pow2nab = 1.0;
   auto strlng = details::stirling();
   for (std::size_t i = 1; i < 10; ++i)
   {
      powNa *= n+a;
      powNb *= n+b;
      pow2nab *= 2*n+a+b;
      sna += strlng[i]/powNa;
      snb += strlng[i]/powNb;
      s2nab += strlng[i]/pow2nab;
   }
   auto sR = sna*snb / s2nab;

   auto t_2 = 0.5*t;
   auto lhs = std::pow(std::sin(t_2), a+0.5)*std::pow(std::cos(t_2), b+0.5);

   // most of the difference from chebfun comes from
   // sumF due to the cos(phi) evaluation
   auto C = std::exp(details::expArgInt(n, a, b))
      * sR*2.0*std::sqrt(2/Math::PI * (n+a)*(n+b)/(2*n+a+b))
      / ((2*n+a+b+1) * lhs );
   auto P = C * sumF;
   auto P_m1 = C * (a+b+2*n)*(a+b+1+2*n)/(4*(a+n)*(b+n)) * sumF_m1;
   // first derivative
   auto cs = JacobiBase::dPnabPnab<Quadrature::natural_t>(n, a, b);
   auto dP = ((cs(1)+cs(2)*t)*P + cs(3)*P_m1)/(-cs(0)*std::sin(t));

   return std::array<MHDFloat, 2>{P, dP};
}

} // namespace JacobiAsy
} // namespace Jacobi
} // namespace Polynomial
} // namespace QuICC
