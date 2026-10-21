/**
 * @file InsulatingTopPerfcondBottomTorShell.cpp
 * @brief Source of the implementation of the mixed boundary insulating top / perfectly conducting bottom in spherical shell
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/InsulatingTopPerfcondBottomTorShell.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

InsulatingTopPerfcondBottomTorShell::InsulatingTopPerfcondBottomTorShell(const Scalar_t lower, const Scalar_t upper,
   const Position pos) :
    ICondition(lower, upper, pos)
{}

InsulatingTopPerfcondBottomTorShell::ACoeff_t InsulatingTopPerfcondBottomTorShell::compute(const int maxN)
{
   ACoeff_t val = ACoeff_t::Zero(maxN + 1);

   if (this->position() == Position::TOP)
   {
      Value bcTop(this->lower(), this->upper(), this->position());
      val = bcTop.compute(maxN);
   }
   else
   {
      const auto ab_a = (-this->a() + this->b()) / this->a() ;
      for (int i = 0; i < val.size(); i++)
      {
         const auto n = static_cast<Scalar_t>(i);
         const auto n2 = n * n;
         val(i) = - ab_a * n2 + 1.0;
         if (i % 2 == 0)
         {
            val(i) *= -1.0;
         }
      }
      // Normalization
      val(0) /= this->c();
   }


   return val;
}

} // namespace Boundary
} // namespace LinearMap
} // namespace Chebyshev
} // namespace SparseSM
} // namespace QuICC
