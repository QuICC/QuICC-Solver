/**
 * @file IProjCrossOperator.cpp
 * @brief Source of the implementation of the generic projection of cross product A ^ B 
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <wigxjpf.h>

// Project includes
//
#include "Types/Math.hpp"
#include "DenseSM/Chebyshev/LinearMap/IProjCrossOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

IProjCrossOperator::IProjCrossOperator(const int rows, const int cols, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB, const Scalar_t lower,
   const Scalar_t upper)
: ILinearMapOperator(rows, cols, lower, upper), mLout(lOut), mMout(mOut), mLa(lA), mMa(mA), mLb(lB), mMb(mB)
{}

MHDFloat IProjCrossOperator::gaunt(const int lA, const int mA, const int lB, const int mB, const int lG, const int mG) const
{
   int lmax = std::max(std::max(lA, lB), lG);

   double val3jA;
   double val3jB;

   wig_table_init(2*lmax,3);
   wig_temp_init(2*lmax);

   /* Note that the arguments to wig3jj, wig6jj and wig9jj are 2*j
    * and 2*m.  To be able to handle half-integer arguments.
    */

   val3jA = wig3jj(2* lA , 2* lB , 2* lG ,
         2* mA, 2* mB , 2* mG);

   val3jB = wig3jj(2* lA , 2* lB , 2* lG ,
         0,  0, 0);

   wig_temp_free();
   wig_table_free();

   MHDFloat Kabg = (std::sqrt(static_cast<MHDFloat>((2*lA + 1)*(2*lB + 1)*(2*lG + 1))))*val3jA*val3jB;

   return Kabg;
}

MHDFloat IProjCrossOperator::elsasser(const int lA, const int mA, const int lB, const int mB, const int lG, const int mG) const
{
   int lmax = std::max(std::max(lA, lB), lG);

   double val3jA;
   double val3jB;

   wig_table_init(2*lmax,3);
   wig_temp_init(2*lmax);

   /* Note that the arguments to wig3jj, wig6jj and wig9jj are 2*j
    * and 2*m.  To be able to handle half-integer arguments.
    */

   val3jA = wig3jj(2* lA , 2* lB , 2* lG ,
         2* mA, 2* mB , 2* mG);

   val3jB = wig3jj(2* lA , 2* (lB + 1) , 2* lG ,
         0,  0, 0);

   wig_temp_free();
   wig_table_free();

   MHDFloat Labg = -(std::sqrt(static_cast<MHDFloat>((2*lA + 1)*(2*lB + 1)*(2*lG + 1)))/2.0)*val3jA*val3jB*std::sqrt(static_cast<MHDFloat>(lA + lB + lG + 2)*static_cast<MHDFloat>(lA + lB - lG + 1)*static_cast<MHDFloat>(lB + lG - lA + 1)*static_cast<MHDFloat>(lA + lG - lB));

   return Labg;
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
