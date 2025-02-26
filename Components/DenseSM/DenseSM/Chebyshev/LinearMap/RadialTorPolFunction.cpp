/**
 * @file RadialTorPolFunction.cpp
 * @brief Source of the implementation of generic radial toroidal/poloidal function
 */

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "Types/Internal/BasicTypes.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

Array RadialTorPolFunction::evaluateDiff(const int p, const Internal::Array& r, const int l,
      const int m, const Internal::MHDFloat lb, const Internal::MHDFloat ub) const
{
   // Evaluate function f on grid
   Matrix f = this->evaluate(r, l, m).cast<MHDFloat>();
   Matrix sf = Utils::computeExpansion(f, this->nN(), lb, ub);
   switch(p)
   {
      case 1:
         f = Utils::evaluateD<1>(sf, this->nN(), lb, ub);
         break;
      case 2:
         f = Utils::evaluateD<2>(sf, this->nN(), lb, ub);
         break;
      case 3:
         f = Utils::evaluateD<3>(sf, this->nN(), lb, ub);
         break;
      case 4:
         f = Utils::evaluateD<4>(sf, this->nN(), lb, ub);
         break;
      default:
         throw std::logic_error("Only derivative order 1-4 are implemented");
         break;
   }

   return f;
}

void RadialTorPolFunction::initParams(const QuICC::Equations::EquationParameters& eqParams)
{
   // nothing here
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
