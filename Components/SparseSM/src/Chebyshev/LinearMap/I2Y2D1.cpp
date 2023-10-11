/**
 * @file I2Y2D1.cpp
 * @brief Source of the implementation of the I^2 Y^2 D sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2D1.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I2Y2D1::I2Y2D1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I2Y2D1::ACoeff_t I2Y2D1::d_3(const ACoeff_t& n) const
   {
      const auto c = Internal::Math::pow(this->a()/2.0, 3);
      return c*(n - 3.0)/(n*(n - 1.0));
   }

   I2Y2D1::ACoeff_t I2Y2D1::d_2(const ACoeff_t& n) const
   {
      const auto a2 = this->a()*this->a();
      const auto& b1 = this->b();
      return a2*b1*(n - 2.0)/(2.0*n*(n - 1.0));
   }

   I2Y2D1::ACoeff_t I2Y2D1::d_1(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto b2 = this->b()*this->b();
      return a1*(a2*n + 3.0*a2 + 4.0*b2*n + 4.0*b2)/(8.0*n*(n + 1.0));
   }

   I2Y2D1::ACoeff_t I2Y2D1::d0(const ACoeff_t& n) const
   {
      const auto& b1 = this->b();
      const auto a2 = this->a()*this->a();
      return a2*b1/((n - 1.0)*(n + 1.0));
   }

   I2Y2D1::ACoeff_t I2Y2D1::d1(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto b2 = this->b()*this->b();
      return -a1*(a2*n - 3.0*a2 + 4.0*b2*n - 4.0*b2)/(8.0*n*(n - 1.0));
   }

   I2Y2D1::ACoeff_t I2Y2D1::d2(const ACoeff_t& n) const
   {
      const auto& b1 = this->b();
      const auto a2 = this->a()*this->a();
      return -a2*b1*(n + 2.0)/(2.0*n*(n + 1.0));
   }

   I2Y2D1::ACoeff_t I2Y2D1::d3(const ACoeff_t& n) const
   {
      const auto c = Internal::Math::pow(this->a()/2.0, 3);
      return -c*(n + 3.0)/(n*(n + 1.0));
   }

   void I2Y2D1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-2, 2, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(7*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -3, ni, this->d_3(n));
         this->convertToTriplets(list, -2, ni, this->d_2(n));
         this->convertToTriplets(list, -1, ni, this->d_1(n));
         this->convertToTriplets(list, 0, ni, this->d0(n));
         this->convertToTriplets(list, 1, ni, this->d1(n));
         this->convertToTriplets(list, 2, ni, this->d2(n));
         this->convertToTriplets(list, 3, ni, this->d3(n));
      }
   }

}
}
}
}
