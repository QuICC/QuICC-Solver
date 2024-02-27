/**
 * @file Utils.cpp
 * @brief Simple utils
 */

// System includes
//
#include "Types/Internal/Typedefs.hpp"
#include <boost/math/special_functions/bessel.hpp>

// Project includes
//
#include "QuICC/SparseSM/Bessel/Value/Utils.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Value {

   void getRoots(std::vector<Internal::MHDFloat>& roots, const int l, const int nRoots)
   {
      using namespace Internal::Literals;
      Internal::MHDFloat nu = static_cast<Internal::MHDFloat>(l) + 0.5_mp;
      boost::math::cyl_bessel_j_zero(nu, 1, nRoots, std::back_inserter(roots));
   }

}
}
}
}
