/**
 * @file InsulatingTopPerfcondBottomPolShell.cpp
 * @brief Source of the implementation of the mixed boundary insulating top / perfectly conducting bottom in spherical shell
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/InsulatingTopPerfcondBottomPolShell.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

InsulatingTopPerfcondBottomPolShell::InsulatingTopPerfcondBottomPolShell(const Scalar_t lower, const Scalar_t upper,
   const Position pos, const int l) :
    ICondition(lower, upper, pos), mL(l)
{}

InsulatingTopPerfcondBottomPolShell::ACoeff_t InsulatingTopPerfcondBottomPolShell::compute(const int maxN)
{
   ACoeff_t val = ACoeff_t::Zero(maxN + 1);

   const auto a_1 = 2.0 / this->a();
   if (this->position() == Position::TOP)
   {
      const auto ab_1 = (this->mL + 1.0) / (this->a() + this->b());
      for (int i = 0; i < val.size(); i++)
      {
         const auto n = static_cast<Scalar_t>(i);
         const auto n2 = n * n;
         val(i) = a_1 * n2 + (ab_1) * this->c();
      }
      // Normalization
      val(0) /= this->c();
   }
   else
   {
      /*
      const auto ab_a2 = (-this->a() + this->b()) / this->a() / this->a() ;
      for (int i = 0; i < val.size(); i++)
      {
         const auto n = static_cast<Scalar_t>(i);
         const auto n2 = n * n;
         const auto n3 = n * n * n;
         val(i) = - a_1 * n2 + ab_a2 * (n2*n2 - n3)/3;
         if (i % 2 == 0)
         {
            val(i) *= -1.0;
         }
      }
      */
      Value bcBot(this->lower(), this->upper(), this->position());
      val = bcBot.compute(maxN);
   }

   return val;
}

} // namespace Boundary
} // namespace LinearMap
} // namespace Chebyshev
} // namespace SparseSM
} // namespace QuICC
