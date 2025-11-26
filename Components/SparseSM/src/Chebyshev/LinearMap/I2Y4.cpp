/**
 * @file I2Y4.cpp
 * @brief Source of the implementation of the I^2 Y^4 sparse operator, with y =
 * ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y4.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

using namespace Internal::Literals;

I2Y4::I2Y4(const int rows, const int cols, const Scalar_t lower,
   const Scalar_t upper) :
    ILinearMapOperator(rows, cols, lower, upper)
{}

I2Y4::ACoeff_t I2Y4::d_6(const ACoeff_t& n) const
{
   return a<6>() / (64.0_mp*(n * (n - 1.0_mp)));
}

I2Y4::ACoeff_t I2Y4::d_5(const ACoeff_t& n) const
{
   return a<5>() * b<1>() / (256.0_mp * n * (n - 1.0_mp));
}

I2Y4::ACoeff_t I2Y4::d_4(const ACoeff_t& n) const
{
   return a<4>() *
          (a<2>() * (n + 2.0_mp) + 12.0_mp * b<2>() * (n + 1.0_mp)) /
          (32.0_mp * n * (n - 1.0_mp) * (n + 1.0_mp));
}

I2Y4::ACoeff_t I2Y4::d_3(const ACoeff_t& n) const
{
   return a<3>() * b<1>() *
          (a<2>() * (n + 3.0_mp) + 4.0_mp * b<2>() * (n + 1.0_mp)) / (8.0_mp * n * (n - 1.0_mp) * (n + 1.0_mp));
}

I2Y4::ACoeff_t I2Y4::d_2(const ACoeff_t& n) const
{
   return -a<2>() * (a<4>() * (n - 5.0_mp) - 48.0_mp * a<2>() * b<2>() - 16.0_mp * b<4>() * (n + 1.0_mp)) /
          (64.0_mp * n * (n - 1.0_mp) * (n + 1.0_mp));
}

I2Y4::ACoeff_t I2Y4::d_1(const ACoeff_t& n) const
{
   return -a<3>() * this->b() * (a<2>() + 2.0_mp * b<2>()) / (4.0_mp * n * (n + 1.0_mp));
}

I2Y4::ACoeff_t I2Y4::d0(const ACoeff_t& n) const
{
   return -a<2>() * (a<4>() + 12.0_mp * a<2>() * b<2>() + 8.0_mp * b<4>()) /
          (16.0_mp * (n - 1.0_mp) * (n + 1.0_mp));
}

I2Y4::ACoeff_t I2Y4::d1(const ACoeff_t& n) const
{
   return this->d_1(n - 1.0_mp);
}

I2Y4::ACoeff_t I2Y4::d2(const ACoeff_t& n) const
{
   return -a<2>() * (a<4>() * (n + 5.0_mp) + 48.0_mp * a<2>() * b<2>() - 16.0_mp * b<4>() * (n - 1.0_mp)) /
          (64.0_mp * n * (n - 1.0_mp) * (n + 1.0_mp));
}

I2Y4::ACoeff_t I2Y4::d3(const ACoeff_t& n) const
{
   return a<3>() * b<1>() * (a<2>() * (n - 3.0_mp) + 4.0_mp * b<2>() * (n - 1.0_mp)) /
          (8.0_mp * n * (n - 1.0_mp) * (n + 1.0_mp));
}

I2Y4::ACoeff_t I2Y4::d4(const ACoeff_t& n) const
{
   return Internal::Math::pow(a<1>() / 2.0_mp, 4) *
          (a<2>() * (n - 2.0_mp) + 12.0_mp * b<2>() * (n - 1.0_mp)) /
          (2.0_mp * n * (n - 1.0_mp) * (n + 1.0_mp));
}

I2Y4::ACoeff_t I2Y4::d5(const ACoeff_t& n) const
{
   return d_5(n + 1.0_mp);
}

I2Y4::ACoeff_t I2Y4::d6(const ACoeff_t& n) const
{
   return d_6(n + 1.0_mp);
}

void I2Y4::buildTriplets(TripletList_t& list) const
{
   ACoeffI ni = ACoeffI::LinSpaced(this->rows() - 2, 2, this->rows() - 1);
   ACoeff_t n = ni.cast<Scalar_t>();

   if (n.size() > 0)
   {
      list.reserve(13 * std::max(this->rows(), this->cols()));
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
   }
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace SparseSM
} // namespace QuICC
