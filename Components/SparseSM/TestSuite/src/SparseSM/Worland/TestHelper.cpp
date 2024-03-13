/**
 * @file TestHelper.cpp
 * @brief Source of test helper
 */

// System includes
//

// Project includes
//
#include "QuICC/TestSuite/SparseSM/Worland/TestHelper.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

namespace Worland {

void setJacobiParameters(Internal::MHDFloat& alpha, Internal::MHDFloat& dBeta,
   const std::string& type)
{
   using namespace QuICC::Internal::Literals;
   if (type == std::string("Chebyshev"))
   {
      alpha = -0.5_mp;
      dBeta = -0.5_mp;
   }
   else if (type == std::string("Legendre"))
   {
      alpha = 0_mp;
      dBeta = -0.5_mp;
   }
   else if (type == std::string("CylEnergy"))
   {
      alpha = 0_mp;
      dBeta = 0_mp;
   }
   else if (type == std::string("SphEnergy"))
   {
      alpha = 0_mp;
      dBeta = 0.5_mp;
   }
   else
   {
      throw std::logic_error("Unknown Worland type");
   }
}

} // namespace Worland
} // namespace SparseSM
} // namespace TestSuite
} // namespace QuICC
