/**
 * @file I3Y4.cpp
 * @brief Source of the implementation of the I^3 Y^4 sparse operator, with y =
 * ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I3Y4.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

using namespace Internal::Literals;

I3Y4::I3Y4(const int rows, const int cols, const Scalar_t lower,
   const Scalar_t upper) :
    ILinearMapOperator(rows, cols, lower, upper)
{}

I3Y4::ACoeff_t I3Y4::d_7(const ACoeff_t& n) const
{
   return a<7>() / (128.0_mp * n * (n - 2.0_mp) * (n - 1.0_mp));
}

I3Y4::ACoeff_t I3Y4::d_6(const ACoeff_t& n) const
{
   return a<6>() * b<1>() / (16.0_mp * n * (n - 2.0_mp) * (n - 1.0_mp));
}

I3Y4::ACoeff_t I3Y4::d_5(const ACoeff_t& n) const
{
   return a<5>() *
          (a<2>() * n + 7.0_mp * a<2>() + 24.0_mp * b<2>() * n +
             24.0_mp * b<2>()) /
          (128.0_mp * n * (n - 2.0_mp) * (n - 1.0_mp) * (n + 1.0_mp));
}

I3Y4::ACoeff_t I3Y4::d_4(const ACoeff_t& n) const
{
   return a<4>() * b<1>() *
          (3.0_mp * a<2>() + 2.0_mp * b<2>() * n + 2.0_mp * b<2>()) /
          (8.0_mp * n * (n - 2.0_mp) * (n - 1.0_mp) * (n + 1.0_mp));
}

I3Y4::ACoeff_t I3Y4::d_3(const ACoeff_t& n) const
{
   return -a<3>() *
          (3.0_mp * a<4>() * n * n - 3.0_mp * a<4>() * n - 30.0_mp * a<4>() +
             24.0_mp * a<2>() * b<2>() * n * n - 72.0_mp * a<2>() * b<2>() * n -
             240.0_mp * a<2>() * b<2>() - 16.0_mp * b<4>() * n * n -
             48.0_mp * b<4>() * n - 32.0_mp * b<4>()) /
          (128.0_mp * n * (n - 2.0_mp) * (n - 1.0_mp) * (n + 1.0_mp) *
             (n + 2.0_mp));
}

I3Y4::ACoeff_t I3Y4::d_2(const ACoeff_t& n) const
{
   return -a<4>() * b<1>() *
          (3.0_mp * a<2>() * n + 9.0_mp * a<2>() + 8.0_mp * b<2>() * n +
             16.0_mp * b<2>()) /
          (16.0_mp * n * (n - 1.0_mp) * (n + 1.0_mp) * (n + 2.0_mp));
}

I3Y4::ACoeff_t I3Y4::d_1(const ACoeff_t& n) const
{
   return -3.0_mp * a<3>() *
          (a<4>() * n + 6.0_mp * a<4>() + 16.0_mp * a<2>() * b<2>() * n +
             64.0_mp * a<2>() * b<2>() + 16.0_mp * b<4>() * n +
             32.0_mp * b<4>()) /
          (128.0_mp * n * (n - 2.0_mp) * (n + 1.0_mp) * (n + 2.0_mp));
}

I3Y4::ACoeff_t I3Y4::d0(const ACoeff_t& n) const
{
   return -3.0_mp * a<4>() * b<1>() * (a<2>() + 2.0_mp * b<2>()) /
          (4.0_mp * (n - 2.0_mp) * (n - 1.0_mp) * (n + 1.0_mp) * (n + 2.0_mp));
}

I3Y4::ACoeff_t I3Y4::d1(const ACoeff_t& n) const
{
   return 3.0_mp * a<3>() *
          (a<4>() * n - 6.0_mp * a<4>() + 16.0_mp * a<2>() * b<2>() * n -
             64.0_mp * a<2>() * b<2>() + 16.0_mp * b<4>() * n -
             32.0_mp * b<4>()) /
          (128.0_mp * n * (n - 2.0_mp) * (n - 1.0_mp) * (n + 2.0_mp));
}

I3Y4::ACoeff_t I3Y4::d2(const ACoeff_t& n) const
{
   return a<4>() * b<1>() *
          (3.0_mp * a<2>() * n - 9.0_mp * a<2>() + 8.0_mp * b<2>() * n -
             16.0_mp * b<2>()) /
          (16.0_mp * n * (n - 2.0_mp) * (n - 1.0_mp) * (n + 1.0_mp));
}

I3Y4::ACoeff_t I3Y4::d3(const ACoeff_t& n) const
{
   return a<3>() *
          (3.0_mp * a<4>() * n * n + 3.0_mp * a<4>() * n - 30.0_mp * a<4>() +
             24.0_mp * a<2>() * b<2>() * n * n + 72.0_mp * a<2>() * b<2>() * n -
             240.0_mp * a<2>() * b<2>() - 16.0_mp * b<4>() * n * n +
             48.0_mp * b<4>() * n - 32.0_mp * b<4>()) /
          (128.0_mp * n * (n - 2.0_mp) * (n - 1.0_mp) * (n + 1.0_mp) *
             (n + 2.0_mp));
}

I3Y4::ACoeff_t I3Y4::d4(const ACoeff_t& n) const
{
   return a<4>() * b<1>() *
          (3.0_mp * a<2>() - 2.0_mp * b<2>() * n + 2.0_mp * b<2>()) /
          (8.0_mp * n * (n - 1.0_mp) * (n + 1.0_mp) * (n + 2.0_mp));
}

I3Y4::ACoeff_t I3Y4::d5(const ACoeff_t& n) const
{
   return -a<5>() *
          (a<2>() * n - 7.0_mp * a<2>() + 24.0_mp * b<2>() * n -
             24.0_mp * b<2>()) /
          (128.0_mp * n * (n - 1.0_mp) * (n + 1.0_mp) * (n + 2.0_mp));
}

I3Y4::ACoeff_t I3Y4::d6(const ACoeff_t& n) const
{
   return -a<6>() * b<1>() / (16.0_mp * n * (n + 1.0_mp) * (n + 2.0_mp));
}

I3Y4::ACoeff_t I3Y4::d7(const ACoeff_t& n) const
{
   return -a<7>() / (128.0_mp * n * (n + 1.0_mp) * (n + 2.0_mp));
}

void I3Y4::buildTriplets(TripletList_t& list) const
{
   ACoeffI ni = ACoeffI::LinSpaced(this->rows() - 3, 3, this->rows() - 1);
   ACoeff_t n = ni.cast<Scalar_t>();

   if (n.size() > 0)
   {
      list.reserve(15 * std::max(this->rows(), this->cols()));
      this->convertToTriplets(list, -7, ni, this->d_7(n));
      this->convertToTriplets(list, -6, ni, this->d_6(n));
      this->convertToTriplets(list, -5, ni, this->d_5(n));
      this->convertToTriplets(list, -4, ni, this->d_4(n));
      this->convertToTriplets(list, -3, ni, this->d_3(n));
      this->convertToTriplets(list, -2, ni, this->d_2(n));
      this->convertToTriplets(list, -1, ni, this->d_1(n));
      this->convertToTriplets(list, 0, ni, this->d0(n));
      this->convertToTriplets(list, 1, ni, this->d1(n));
      this->convertToTriplets(list, 2, ni, this->d2(n));
      this->convertToTriplets(list, 3, ni, this->d3(n));
      this->convertToTriplets(list, 4, ni, this->d4(n));
      this->convertToTriplets(list, 5, ni, this->d5(n));
      this->convertToTriplets(list, 6, ni, this->d6(n));
      this->convertToTriplets(list, 7, ni, this->d7(n));
   }
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace SparseSM
} // namespace QuICC
