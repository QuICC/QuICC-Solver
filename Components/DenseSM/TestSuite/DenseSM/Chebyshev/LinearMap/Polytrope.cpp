/**
 * @file Polytrope.cpp
 * @brief Source of the implementation of the polytropic density field
 */

// System includes
//
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <iostream>
// Project includes
//
#include "TestSuite/DenseSM/Chebyshev/LinearMap/Polytrope.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

Polytrope::Polytrope(const MHDFloat rratio, const MHDFloat Nrho, const MHDFloat npoly) :
   mRratio(rratio), mNrho(Nrho), mNpoly(npoly)
{}

int Polytrope::nN() const
{
   return 2;
}

Internal::Array Polytrope::evaluate(const Internal::Array& r, const int l, const int m) const
{
   using namespace Internal::Literals;

   Internal::Array val;
   if(l != 0)
   {
      val = 0*r;

   }
   else
   {
      
      const Internal::MHDFloat zetaO = (mRratio+1)/( mRratio * Internal::Math::exp(mNrho/mNpoly) + 1);

      const Internal::MHDFloat c0 = (2.0_mp*zetaO - mRratio - 1.0_mp ) / (1.0_mp - mRratio);

      const Internal::MHDFloat c1 = (1.0_mp + mRratio)*(1.0_mp - zetaO) / (1.0_mp - mRratio) / (1.0_mp - mRratio);

      const Internal::MHDFloat zetaI = (1.0_mp + mRratio - zetaO)/mRratio;

      auto zeta = c0 + c1/r.array();

      val = zeta.array().pow(mNpoly);

      //val = (5.602043293950358 * (1.7500004377835053*1e-7 +  0.4224998862499716 * r.array()).pow(2))/r.array().pow(2);
      
      //Another test
      //val = 5.6*(1*1e-7 +  r.array()).pow(2) / r.array().pow(2);

      
      /*
      //dipolar, to test:
      Internal::MHDFloat c1 = (10.0_mp/13.0_mp)*Internal::Math::sqrt(15.0_mp/732251.0_mp);
      val = c1*(49.0_mp - 566.0_mp*r.array() + 400.0_mp*r.array().pow(2));
      */

      //std::cerr << "rad = "<< r <<  " \n"; //ok
      //std::cerr << "density = "<< val <<  " \n"; //ok

      //throw std::logic_error("l=0 case");
   }

   return val;
}

std::vector<int> Polytrope::ls() const
{
   std::vector<int> l = {1};

   return l;
}

}
}
}
}
} // namespace QuICC
