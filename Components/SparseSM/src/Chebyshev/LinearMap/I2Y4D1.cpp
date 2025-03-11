/**
 * @file I2Y4D1.cpp
 * @brief Source of the implementation of the I^2 Y^4 D sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y4D1.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I2Y4D1::I2Y4D1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I2Y4D1::ACoeff_t I2Y4D1::d_5(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a5 = a2*a2*a1;

      return a5*(n - 5.0)/(32.0*n*(n - 1.0));
   }


   I2Y4D1::ACoeff_t I2Y4D1::d_4(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      return a4*b1*(n - 4.0)/(4.0*n*(n - 1.0));
   }

   I2Y4D1::ACoeff_t I2Y4D1::d_3(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      return a3*(n - 3.0)*(3.0*a2*n + 5*a2 + 24*b2*(n + 1))/(32.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y4D1::ACoeff_t I2Y4D1::d_2(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      return a2*b1*(n - 2.0)*(a2*(n + 2.0) + 2.0*b2*(n + 1.0))/(2.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y4D1::ACoeff_t I2Y4D1::d_1(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto b4 = b2*b2;
      return a1*(a4*(n + 5.0) + 12.0*a2*b2*(n + 3.0) + 8.0*b4*(n + 1.0))/(16.0*n*(n + 1.0));
   }

   I2Y4D1::ACoeff_t I2Y4D1::d0(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      return a2*b1*(3.0*a2 + 4.0*b2)/(2.0*(n - 1.0)*(n + 1.0));
   }

   I2Y4D1::ACoeff_t I2Y4D1::d1(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto b4 = b2*b2;
      return -a1*(a4*(n - 5.0) + 12.0*a2*b2*(n - 3.0) + 8*b4*(n - 1.0))/(16.0*n*(n - 1.0));
   }

   I2Y4D1::ACoeff_t I2Y4D1::d2(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      return -a2*b1*(n + 2.0)*(a2*(n - 2.0) + 2.0*b2*(n - 1.0))/(2.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y4D1::ACoeff_t I2Y4D1::d3(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto b2 = this->b()*this->b();
      return -a3*(n + 3.0)*(3.0*a2*n - 5.0*a2 + 24*b2*(n - 1.0))/(32.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y4D1::ACoeff_t I2Y4D1::d4(const ACoeff_t& n) const
   {
      const auto& b1 = this->b();
      const auto a2 = this->a()*this->a();
      const auto a4 = a2*a2;
      return -a4*b1*(n + 4.0)/(4.0*n*(n + 1.0));
   }

   I2Y4D1::ACoeff_t I2Y4D1::d5(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a5 = a1*a2*a2;
      return -a5*(n + 5.0)/(32.0*n*(n + 1.0));
   }

   void I2Y4D1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-2, 2, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(11*std::max(this->rows(),this->cols()));
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
      }
   }

}
}
}
}
