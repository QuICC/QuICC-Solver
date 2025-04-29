/**
 * @file I2Y2D1Y1.cpp
 * @brief Source of the implementation of the I^2 Y D Y sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2D1Y1.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I2Y2D1Y1::I2Y2D1Y1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I2Y2D1Y1::~I2Y2D1Y1()
   {
   }

   I2Y2D1Y1::ACoeff_t I2Y2D1Y1::d_4(const ACoeff_t& n) const
   {
      return Internal::Math::pow(this->a(),4)*(n - 3.0)/(16*n*(n - 1.0));
   }

   I2Y2D1Y1::ACoeff_t I2Y2D1Y1::d_3(const ACoeff_t& n) const
   {
      return Internal::Math::pow(this->a(),3)*(this->b())*(3*n - 7.0)/(8*n*(n - 1.0));
   }

   I2Y2D1Y1::ACoeff_t I2Y2D1Y1::d_2(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return a2*(a2*n*n - 3*a2 + 6*b2*n*n - 4*b2*n - 10*b2)/(8*n*(n - 1)*(n + 1));
   }

   I2Y2D1Y1::ACoeff_t I2Y2D1Y1::d_1(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return (this->a())*(this->b())*(3*a2*n + 7*a2 + 4*b2*n + 4*b2)/(8*n*(n + 1));
   }

   I2Y2D1Y1::ACoeff_t I2Y2D1Y1::d0(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return a2*(a2 + 4*b2)/(4*(n - 1)*(n + 1));
   }

   I2Y2D1Y1::ACoeff_t I2Y2D1Y1::d1(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -(this->a())*(this->b())*(3*a2*n - 7*a2 + 4*b2*n - 4*b2)/(8*n*(n - 1));
   }

   I2Y2D1Y1::ACoeff_t I2Y2D1Y1::d2(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -a2*(a2*n*n - 3*a2 + 6*b2*n*n + 4*b2*n - 10*b2)/(8*n*(n - 1)*(n + 1));
   }

   I2Y2D1Y1::ACoeff_t I2Y2D1Y1::d3(const ACoeff_t& n) const
   {
      auto a3 = this->a()*this->a()*this->a();
      return -a3*(this->b())*(3*n + 7)/(8*n*(n + 1));
   }

   I2Y2D1Y1::ACoeff_t I2Y2D1Y1::d4(const ACoeff_t& n) const
   {
      auto a4 = this->a()*this->a()*this->a()*this->a();
      return -a4*(n + 3)/(16*n*(n + 1));
   }

   void I2Y2D1Y1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-2, 2, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(9*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -4, ni, this->d_4(n));
         this->convertToTriplets(list, -3, ni, this->d_3(n));
         this->convertToTriplets(list, -2, ni, this->d_2(n));
         this->convertToTriplets(list, -1, ni, this->d_1(n));
         this->convertToTriplets(list, 0, ni, this->d0(n));
         this->convertToTriplets(list, 1, ni, this->d1(n));
         this->convertToTriplets(list, 2, ni, this->d2(n));
         this->convertToTriplets(list, 3, ni, this->d3(n));
         this->convertToTriplets(list, 4, ni, this->d4(n));
      }
   }

}
}
}
}
