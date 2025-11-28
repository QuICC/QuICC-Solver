/**
 * @file ShellPolytrope.cpp
 * @brief Source of the implementation of the polytropic density field
 */

// System includes
//

// Project includes
//
#include "QuICC/TestSuite/Framework/Io/ShellPolytrope.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Io {

ShellPolytrope::ShellPolytrope(const MHDFloat rratio, const MHDFloat Nrho,
   const MHDFloat npoly) :
    mRratio(rratio), mNrho(Nrho), mNpoly(npoly)
{}

int ShellPolytrope::nN() const
{
   return 2;
}

Internal::Array ShellPolytrope::evaluate(const Internal::Array& r, const int l,
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
      const Internal::MHDFloat rratio = static_cast<Internal::MHDFloat>(mRratio);
      const Internal::MHDFloat nRho = static_cast<Internal::MHDFloat>(mNrho);
      const Internal::MHDFloat nPoly = static_cast<Internal::MHDFloat>(mNpoly);

      const Internal::MHDFloat zetaO =
         (rratio + 1) / (rratio * Internal::Math::exp(nRho / nPoly) + 1);

      const Internal::MHDFloat c0 =
         (2.0_mp * zetaO - rratio - 1.0_mp) / (1.0_mp - rratio);

      const Internal::MHDFloat c1 = (1.0_mp + rratio) * (1.0_mp - zetaO) /
                                    (1.0_mp - rratio) / (1.0_mp - rratio);

      //const Internal::MHDFloat zetaI = (1.0_mp + rratio - zetaO) / rratio;

      auto zeta = c0 + c1 / r.array();

      val = zeta.array().pow(nPoly);
   }

   return val;
}

std::vector<int> ShellPolytrope::ls() const
{
   std::vector<int> l = {1};

   return l;
}

} // namespace Io
} // namespace Framework
} // namespace TestSuite
} // namespace QuICC
