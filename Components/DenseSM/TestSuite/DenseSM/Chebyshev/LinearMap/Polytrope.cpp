/**
 * @file Polytrope.cpp
 * @brief Source of the implementation of the polytropic density field
 */

// System includes
//
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <sstream>
// Project includes
//
#include "TestSuite/DenseSM/Chebyshev/LinearMap/Polytrope.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

Polytrope::Polytrope(const MHDFloat rratio, const MHDFloat Nrho,
   const MHDFloat npoly) :
    mRratio(rratio), mNrho(Nrho), mNpoly(npoly)
{}

int Polytrope::nN() const
{
   return 2;
}

Internal::Array Polytrope::evaluate(const Internal::Array& r, const int l,
   const int m) const
{
   using namespace Internal::Literals;

   Internal::Array val;
   if (l != 0)
   {
      val = 0 * r;
   }
   else
   {

      const Internal::MHDFloat zetaO =
         (mRratio + 1) / (mRratio * Internal::Math::exp(mNrho / mNpoly) + 1);

      const Internal::MHDFloat c0 =
         (2.0_mp * zetaO - mRratio - 1.0_mp) / (1.0_mp - mRratio);

      const Internal::MHDFloat c1 = (1.0_mp + mRratio) * (1.0_mp - zetaO) /
                                    (1.0_mp - mRratio) / (1.0_mp - mRratio);

      const Internal::MHDFloat zetaI = (1.0_mp + mRratio - zetaO) / mRratio;

      auto zeta = c0 + c1 / r.array();

      val = zeta.array().pow(mNpoly);
   }

   return val;
}

std::vector<int> Polytrope::ls() const
{
   std::vector<int> l = {1};

   return l;
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace TestSuite
} // namespace QuICC
