/**
 * @file Utils.cpp
 * @brief Source of the utils
 * product A ^ B
 */

// System includes
//
#include <cassert>
#include <cmath>
#include <wigxjpf.h>

// Project includes
//
#include "DenseSM/Utils.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace DenseSM {

namespace Utils {

MHDFloat gaunt(const int lA, const int mA, const int lB,
   const int mB, const int lG, const int mG)
{
   int lmax = std::max(std::max(lA, lB), lG);

   double val3jA;
   double val3jB;

   wig_table_init(2 * lmax, 3);
   wig_temp_init(2 * lmax);

   /* Note that the arguments to wig3jj, wig6jj and wig9jj are 2*j
    * and 2*m.  To be able to handle half-integer arguments.
    */

   val3jA = wig3jj(2 * lA, 2 * lB, 2 * lG, 2 * mA, 2 * mB, -2 * mG);

   val3jB = wig3jj(2 * lA, 2 * lB, 2 * lG, 0, 0, 0);

   wig_temp_free();
   wig_table_free();

   MHDFloat Kabg = 0;
   if (val3jA != 0 && val3jB != 0)
   {
      MHDFloat ca = static_cast<MHDFloat>(2 * lA + 1);
      MHDFloat cb = static_cast<MHDFloat>(2 * lB + 1);
      MHDFloat cg = static_cast<MHDFloat>(2 * lG + 1);
      Kabg = std::sqrt(ca * cb * cg / (4.0 * Math::PI)) * val3jA * val3jB;
   }

   return Kabg;
}

MHDFloat elsasser(const int lA, const int mA, const int lB,
   const int mB, const int lG, const int mG)
{
   int lmax = std::max(std::max(lA, lB + 1), lG);

   double val3jA;
   double val3jB;

   wig_table_init(2 * lmax, 3);
   wig_temp_init(2 * lmax);

   /* Note that the arguments to wig3jj, wig6jj and wig9jj are 2*j
    * and 2*m.  To be able to handle half-integer arguments.
    */

   val3jA = wig3jj(2 * lA, 2 * lB, 2 * lG, 2 * mA, 2 * mB, -2 * mG);

   val3jB = wig3jj(2 * lA, 2 * (lB + 1), 2 * lG, 0, 0, 0);

   wig_temp_free();
   wig_table_free();

   MHDFloat Labg = 0;
   if (val3jA != 0 && val3jB != 0)
   {
      MHDFloat ca = static_cast<MHDFloat>(2 * lA + 1);
      MHDFloat cb = static_cast<MHDFloat>(2 * lB + 1);
      MHDFloat cg = static_cast<MHDFloat>(2 * lG + 1);
      MHDFloat clabg2 = static_cast<MHDFloat>(lA + lB + lG + 2);
      MHDFloat clab_g1 = static_cast<MHDFloat>(lA + lB - lG + 1);
      MHDFloat clbg_a1 = static_cast<MHDFloat>(lB + lG - lA + 1);
      MHDFloat clag_b = static_cast<MHDFloat>(lA + lG - lB);
      Labg = (std::sqrt(ca * cb * cg / (4.0 * Math::PI)) / 2.0) * val3jA *
             val3jB * std::sqrt(clabg2 * clab_g1 * clbg_a1 * clag_b);
   }

   return Labg;
}

} // namespace Utils
} // namespace DenseSM
} // namespace QuICC
