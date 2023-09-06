/** 
 * @file I2Y3.cpp
 * @brief Source of the implementation of the I^2 Y^3 sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y3.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I2Y3::I2Y3(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I2Y3::ACoeff_t I2Y3::d_5(const ACoeff_t& n) const
   {
      const auto c = precision::pow(this->a()/2.0,5);
      return c/(n*(n - 1.0));
   }

   I2Y3::ACoeff_t I2Y3::d_4(const ACoeff_t& n) const
   {
      const auto& b1 = this->b(); 
      const auto c = precision::pow(this->a()/2.0,4);
      return c*3.0*b1/(n*(n - 1.0));
   }

   I2Y3::ACoeff_t I2Y3::d_3(const ACoeff_t& n) const
   {
      const auto a2 = this->a()*this->a(); 
      const auto b2 = this->b()*this->b(); 
      const auto c = precision::pow(this->a()/2.0,3);
      return c*(a2*(n + 3.0) + 12.0*b2*(n + 1.0))/(4.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y3::ACoeff_t I2Y3::d_2(const ACoeff_t& n) const
   {
      const auto& b1 = this->b(); 
      const auto a2 = this->a()*this->a();
      const auto b2 = this->b()*this->b();
      const auto c = precision::pow(this->a()/2.0, 2);
      return c*b1*(3.0*a2 + 2.0*b2*(n + 1.0))/(2.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y3::ACoeff_t I2Y3::d_1(const ACoeff_t& n) const
   {
      const auto a2 = this->a()*this->a();
      const auto a3 = this->a()*a2;
      const auto b2 = this->b()*this->b();
      return -a3*(a2 + 6.0*b2)/(16.0*n*(n + 1.0));
   }

   I2Y3::ACoeff_t I2Y3::d0(const ACoeff_t& n) const
   {
      const auto& b1 = this->b();
      const auto a2 = this->a()*this->a();
      const auto b2 = this->b()*this->b();
      return -a2*b1*(3.0*a2 + 4.0*b2)/(8.0*(n - 1.0)*(n + 1.0));
   }

   I2Y3::ACoeff_t I2Y3::d1(const ACoeff_t& n) const
   {
      return this->d_1(n - 1.0);
   }

   I2Y3::ACoeff_t I2Y3::d2(const ACoeff_t& n) const
   {
      const auto& b1 = this->b();
      const auto a2 = this->a()*this->a();
      const auto b2 = this->b()*this->b();
      return -a2*b1*(3.0*a2 - 2.0*b2*(n - 1.0))/(8.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y3::ACoeff_t I2Y3::d3(const ACoeff_t& n) const
   {
      const auto a2 = this->a()*this->a(); 
      const auto b2 = this->b()*this->b(); 
      const auto c = precision::pow(this->a()/2.0,3);
      return c*(a2*(n - 3.0) + 12.0*b2*(n - 1.0))/(4.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y3::ACoeff_t I2Y3::d4(const ACoeff_t& n) const
   {
      const auto& b1 = this->b();
      const auto c = precision::pow(this->a()/2.0,4);
      return c*3.0*b1/(n*(n + 1.0));
   }

   I2Y3::ACoeff_t I2Y3::d5(const ACoeff_t& n) const
   {
      const auto c = precision::pow(this->a()/2.0,5);
      return c/(n*(n + 1.0));
   }

   void I2Y3::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-2, 2, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(9*std::max(this->rows(),this->cols()));
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

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
