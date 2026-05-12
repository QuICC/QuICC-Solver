/**
 * @file OrthogonalValueD2.cpp
 * @brief Source of the implementation of orthogonal Galerkin stencil for value and second derivative boundary condition
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/Stencil/OrthogonalValueD2.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

namespace Stencil {

OrthogonalValueD2::OrthogonalValueD2(const int rows, const int cols,
   const Scalar_t alpha, const Scalar_t dBeta, const int l) :
    IStencilOperator(rows, cols, alpha, dBeta, l)
{
}

void OrthogonalValueD2::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (2 * (this->rows() + 3) + this->mL);

   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);

   Polynomial::Worland::Wnl bW(this->mcAlpha + 3, Polynomial::Worland::worland_default_t::DBETA);
   Internal::Matrix tmp(igrid.size(), this->rows() + 3);
   bW.compute<Internal::MHDFloat>(tmp, tmp.cols(), this->mL, igrid,
      Internal::Array(), ev::Set());
   Internal::Matrix opBwd(igrid.size(), this->rows()-2);
   // Set first basis function
   const Scalar_t& a = this->mcAlpha;
   const Scalar_t& l = this->mL;
   Scalar_t cnst = (99+4*l*l+40*a+4*a*a+8*l*(5+a))/(32.*std::pow(l,0.75));
   opBwd.col(0) = cnst*(1.0 - igrid.array().pow(2)) * igrid.array().pow(l) * (-(5. + 2*l) + (1 + 2*l)*igrid.array().pow(2));
   // Build linear combination
   for(int i = 1; i < opBwd.cols(); i++)
   {
      opBwd.col(i) = d1(i+1, this->mL - 1)*tmp.col(i+2) + d2(i+1, this->mL - 1)*tmp.col(i+1) + d3(i+1, this->mL - 1)*tmp.col(i) + d4(i+1, this->mL - 1)*tmp.col(i-1);
   }

   Polynomial::Worland::Wnl W;
   Internal::Matrix opFwd(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwd, opFwd.cols(), this->mL, igrid,
      iweights, ev::Set());

   mat = opFwd.transpose() * opBwd;
   Internal::Array norm = ((mat.transpose() * mat)).diagonal().array().sqrt().pow(-1).matrix();

   mat = mat * norm.asDiagonal();
}

OrthogonalValueD2::Scalar_t OrthogonalValueD2::c1(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   const Scalar_t a2 = a*a;
   const Scalar_t a3 = a2*a;
   const Scalar_t l2 = l*l;
   const Scalar_t n2 = n*n;
   const Scalar_t n3 = n2*n;
   const Scalar_t n4 = n2*n2;
   Scalar_t c = n*(n+1)*(2*l+9+2*a+2*n)*(2*l+5+2*a+4*n)*(36-200*n+4*l2*a3-96*l+27*a+16*n*l*a3-92*l*a+320*n2-200*n*a+240*n*l+16*n2*a3+8*a2+372*n2*a-68*n*a2+320*n3+320*l*a*n+48*l2+a3+128*n*a2*l+304*n2*a*l+32*n3*l*a+48*n2*l*a2+96*l2*a*n+16*n*l2*a2+16*n2*l2*a+448*n2*l+128*n*l2-32*l*a2+76*l2*a-8*n*a3+128*n3*l+208*n3*a+136*n2*a2-4*a3*l+16*n4*a+64*n2*l2+32*n3*a2+32*l2*a2+64*n4)*(2*l+3+2*a+4*n);

   return c;
}

OrthogonalValueD2::Scalar_t OrthogonalValueD2::c2(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   const Scalar_t a2 = a*a;
   const Scalar_t a3 = a2*a;
   const Scalar_t a4 = a2*a2;
   const Scalar_t a5 = a3*a2;
   const Scalar_t l2 = l*l;
   const Scalar_t l3 = l2*l;
   const Scalar_t n2 = n*n;
   const Scalar_t n3 = n2*n;
   const Scalar_t n4 = n2*n2;
   const Scalar_t n5 = n3*n2;
   const Scalar_t n6 = n3*n3;
   Scalar_t c = -(2*l+9+2*a+4*n)*(-1260+28*l*a4+4764*n+1700*l2*a3+648*l-2883*a+7386*n*l*a3+16*l3*a4+56*n*a4+1746*l*a+31536*n2+6253*n*a+25272*n*l+7854*n2*a3+264*l2*a4+5352*l*a3*n2+4624*n3*l2*a+1120*n2*l3*a+384*n3*l3+3648*n3*a3+1120*n2*a4+16*a5*l2-2006*a2+46966*n2*a+3116*n*a2+38160*n3+40366*l*a*n+96*n3*l3*a+3312*l2+864*l3-605*a3-82*a4+24744*n*a2*l+62400*n2*a*l+31368*n3*l*a+2696*l3*n*a+27672*n2*l*a2+22316*l2*a*n+10576*n*l2*a2+17784*n2*l2*a+51264*n2*l+832*n3*a3*l+256*n3*a4+416*n4*a3+10752*n4*l+16656*n*l2+384*n2*a4*l+64*n2*a5+64*n*a5*l+1188*l*a2+8896*n3*a2*l-4*a5+4528*n4*a2+320*n5*a2+7260*l2*a+640*n3*a2*l2+800*n4*a2*l+504*l2*a3*n2+5232*n2*a2*l2+864*n*a2*l3+160*l2*a4*n+160*n2*a2*l3+687*n*a3+36192*n3*l+44844*n3*a+200*l3*a3+4416*n5+88*l3*n*a3+384*n6+5888*n4*l*a+96*n6*a+1152*n4*l2+2384*n5*a+1152*n5*l+27364*n2*a2+310*a3*l+880*l3*a2+8256*n3*l2+16280*n4*a+288*n5*l*a+288*n4*l2*a+19584*n2*l2+1920*l3*n2+19328*n3*a2+1560*l3*a+5208*l2*a2+2592*n*l3+2148*l2*a3*n+1088*l*a4*n+19200*n4)*n*(2*l+3+2*a+4*n);

   return c;
}

OrthogonalValueD2::Scalar_t OrthogonalValueD2::c3(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   const Scalar_t a2 = a*a;
   const Scalar_t a3 = a2*a;
   const Scalar_t a4 = a2*a2;
   const Scalar_t a5 = a3*a2;
   const Scalar_t l2 = l*l;
   const Scalar_t l3 = l2*l;
   const Scalar_t n2 = n*n;
   const Scalar_t n3 = n2*n;
   const Scalar_t n4 = n2*n2;
   const Scalar_t n5 = n3*n2;
   const Scalar_t n6 = n3*n3;
   Scalar_t c = (n+a+3)*(2*l+11+2*a+4*n)*(126*l*a4+3132*n+904*l2*a3-1620*a+4146*n*l*a3+8*l3*a4+260*n*a4+1800*l*a+13056*n2+6405*n*a+9240*n*l+4398*n2*a3+140*l2*a4+2760*l*a3*n2+3728*n3*l2*a+896*n2*l3*a+384*n3*l3+1888*n3*a3+608*n2*a4+8*a5*l2-1629*a2+23514*n2*a+4728*n*a2+18960*n3+19590*l*a*n+96*n3*l3*a-606*a3-99*a4+13480*n*a2*l+31800*n2*a*l+19912*n3*l*a+1640*l3*n*a+14400*n2*l*a2+10540*l2*a*n+5312*n*l2*a2+11112*n2*l2*a+16*a5*n+8*a5*l+24480*n2*l+512*n3*a3*l+128*n3*a4+256*n4*a3+8832*n4*l+6672*n*l2+192*n2*a4*l+32*n2*a5+32*n*a5*l+1890*l*a2+5568*n3*a2*l-6*a5+2848*n4*a2+256*n5*a2+2640*l2*a+512*n3*a2*l2+640*n4*a2*l+312*l2*a3*n2+3264*n2*a2*l2+544*n*a2*l3+80*l2*a4*n+128*n2*a2*l3+1615*n*a3+23136*n3*l+23276*n3*a+96*l3*a3+3648*n5+56*l3*n*a3+384*n6+4768*n4*l*a+96*n6*a+1152*n4*l2+1936*n5*a+1152*n5*l+14976*n2*a2+736*a3*l+376*l3*a2+6720*n3*l2+10440*n4*a+288*n5*l*a+288*n4*l2*a+12096*n2*l2+1536*l3*n2+10144*n3*a2+480*l3*a+2548*l2*a2+1440*n*l3+1092*l2*a3*n+592*l*a4*n+12480*n4)*(2*l+5+2*a+4*n);

   return c;
}

OrthogonalValueD2::Scalar_t OrthogonalValueD2::c4(const int n, const int l) const
{
   const Scalar_t& a = this->mcAlpha;
   const Scalar_t a2 = a*a;
   const Scalar_t a3 = a2*a;
   const Scalar_t l2 = l*l;
   const Scalar_t n2 = n*n;
   const Scalar_t n3 = n2*n;
   const Scalar_t n4 = n2*n2;
   Scalar_t c = -(2*n-1+2*l)*(n+a+3)*(n+2+a)*(2*l+11+2*a+4*n)*(2*l+9+2*a+4*n)*(540+1656*n+4*l2*a3+720*l+423*a+16*n*l*a3+564*l*a+1664*n2+1232*n*a+1520*n*l+16*n2*a3+108*a2+1092*n2*a+300*n*a2+576*n3+1024*l*a*n+240*l2+9*a3+224*n*a2*l+400*n2*a*l+32*n3*l*a+48*n2*l*a2+128*l2*a*n+16*n*l2*a2+16*n2*l2*a+832*n2*l+256*n*l2+144*l*a2+188*l2*a+24*n*a3+128*n3*l+272*n3*a+232*n2*a2+12*a3*l+16*n4*a+64*n2*l2+32*n3*a2+48*l2*a2+64*n4);

   return c;
}

OrthogonalValueD2::Scalar_t OrthogonalValueD2::d1(const int n_, const int l_) const
{
   if(this->mcAlpha != -0.5)
   {
      throw std::logic_error("Only implemented for alpha = -1/2");
   }

   const Scalar_t l = static_cast<Scalar_t>(l_);
   const Scalar_t n = static_cast<Scalar_t>(n_);
   const Scalar_t l2 = l*l;
   Scalar_t num = 195+28*l2*(1+4*n)*(5+4*n)+16*n*(2+n)*(-29+28*n*(2+n))+4*l*(-115+4*n*(55+14*n*(11+4*n)));
   Scalar_t den = 2*(3+l+n)*(5+2*n)*(1+l+2*n)*std::pow(5+l+2*n,2)*(6+l+2*n);

   return num/den;
}

OrthogonalValueD2::Scalar_t OrthogonalValueD2::d2(const int n_, const int l_) const
{
   if(this->mcAlpha != -0.5)
   {
      throw std::logic_error("Only implemented for alpha = -1/2");
   }

   const Scalar_t l = static_cast<Scalar_t>(l_);
   const Scalar_t n = static_cast<Scalar_t>(n_);
   const Scalar_t l2 = l*l;
   const Scalar_t l3 = l2*l;
   Scalar_t num = -((4+l+2*n)*(-1995+56*l3*(8+3*n)*(1+4*n)*(5+4*n)+l2*(6300+4*n*(15767+2*n*(11937+56*n*(109+18*n))))+n*(18673+2*n*(55921+8*n*(10065+2*n*(3035+14*n*(59+6*n)))))+2*l*(140+n*(41671+4*n*(26337+4*n*(5657+14*n*(143+18*n)))))));
   Scalar_t den = 2*(3+l+n)*(5+2*n)*(1+l+2*n)*(2+l+2*n)*std::pow(5+l+2*n,2)*std::sqrt((1+n)*(4+l+n)*(7+2*n)*(4+l+2*n)*(6+l+2*n)*(3+2*l+2*n));

   return num/den;
}

OrthogonalValueD2::Scalar_t OrthogonalValueD2::d3(const int n_, const int l_) const
{
   if(this->mcAlpha != -0.5)
   {
      throw std::logic_error("Only implemented for alpha = -1/2");
   }

   const Scalar_t l = static_cast<Scalar_t>(l_);
   const Scalar_t n = static_cast<Scalar_t>(n_);
   const Scalar_t l2 = l*l;
   const Scalar_t l3 = l2*l;
   Scalar_t num = 3780+28*l3*(5+4*n)*(9+4*n)*(-1+6*n)+4*l2*(-1575+n*(5197+2*n*(7317+56*n*(89+18*n))))+n*(7403+2*n*(18121+8*n*(4815+2*n*(1985+14*n*(49+6*n)))))+l*(-4095+2*n*(9331+4*n*(11847+4*n*(3627+28*n*(59+9*n)))));
   Scalar_t den = 2*std::pow(1+l+2*n,2)*(5+l+2*n)*std::sqrt(n*(1+n)*(3+l+n)*(4+l+n)*(5+2*n)*(7+2*n)*(2+l+2*n)*(6+l+2*n)*(1+2*l+2*n)*(3+2*l+2*n));

   return num/den;
}

OrthogonalValueD2::Scalar_t OrthogonalValueD2::d4(const int n_, const int l_) const
{
   if(this->mcAlpha != -0.5)
   {
      throw std::logic_error("Only implemented for alpha = -1/2");
   }

   const Scalar_t l = static_cast<Scalar_t>(l_);
   const Scalar_t n = static_cast<Scalar_t>(n_);
   const Scalar_t l2 = l*l;
   Scalar_t num = -((4+l+2*n)*std::sqrt((-1+n)*(2+l+n)*(3+2*n)*(-1+2*l+2*n))*(2835+28*l2*(5+4*n)*(9+4*n)+16*n*(4+n)*(139+28*n*(4+n))+4*l*(945+4*n*(531+14*n*(23+4*n)))));
   Scalar_t den = 2*std::pow(1+l+2*n,2)*(2+l+2*n)*(5+l+2*n)*std::sqrt(n*(1+n)*(3+l+n)*(4+l+n)*(5+2*n)*(7+2*n)*(l+2*n)*(6+l+2*n)*(1+2*l+2*n)*(3+2*l+2*n));

   return num/den;
}

} // namespace Stencil
} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
