/**
 * @file I4Y4D1.cpp
 * @brief Source of the implementation of the I^4 Y^4 D sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4D1.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I4Y4D1::I4Y4D1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I4Y4D1::ACoeff_t I4Y4D1::d_7(const ACoeff_t& n) const
   {
      const auto c = Internal::Math::pow(this->a()/2.0, 7);
      return c*(n - 7.0)/(n*(n - 3.0)*(n - 2.0)*(n - 1.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d_6(const ACoeff_t& n) const
   {
      const auto a2 = this->a()*this->a();
      const auto a6 = a2*a2*a2;
      const auto& b1 = this->b();
      return a6*b1*(n - 6.0)/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d_5(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a5 = a2*a2*a1;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      return a5*(n - 5.0)*(a2*n + 13.0*a2 + 24.0*b2*n + 24.0*b2)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d_4(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      return a4*b1*(n - 4.0)*(3.0*a2 + b2*n + b2)/(4.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d_3(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto b4 = b2*b2;
      return -a3*(3.0*a4*n.pow(2) - 15.0*a4*n - 102.0*a4 + 24.0*a2*b2*n.pow(2) - 216.0*a2*b2*n - 528.0*a2*b2 - 16.0*b4*n.pow(2) - 48.0*b4*n - 32.0*b4)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d_2(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      return -a4*b1*(3.0*a2*n.pow(2) - 3.0*a2*n - 78.0*a2 + 8.0*b2*n.pow(2) - 24.0*b2*n - 80.0*b2)/(16.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d_1(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto b4 = b2*b2;
      return -3.0*a3*(a4*n.pow(3) + 10.0*a4*n.pow(2) - 29.0*a4*n - 190.0*a4 + 16.0*a2*b2*n.pow(3) + 96.0*a2*b2*n.pow(2) - 304.0*a2*b2*n - 1344.0*a2*b2 + 16.0*b4*n.pow(3) + 32.0*b4*n.pow(2) - 144.0*b4*n - 288.0*b4)/(128.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d0(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      return -3.0*a4*b1*(a2*n.pow(2) - 14.0*a2 + 2.0*b2*n.pow(2) - 18.0*b2)/(2*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d1(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto b4 = b2*b2;
      return 3.0*a3*(a4*n.pow(3) - 10.0*a4*n.pow(2) - 29.0*a4*n + 190.0*a4 + 16.0*a2*b2*n.pow(3) - 96.0*a2*b2*n.pow(2) - 304.0*a2*b2*n + 1344.0*a2*b2 + 16.0*b4*n.pow(3) - 32.0*b4*n.pow(2) - 144.0*b4*n + 288.0*b4)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d2(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      return a4*b1*(3.0*a2*n.pow(2) + 3.0*a2*n - 78.0*a2 + 8.0*b2*n.pow(2) + 24.0*b2*n - 80.0*b2)/(16.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d3(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto b4 = b2*b2;
      return a3*(3.0*a4*n.pow(2) + 15.0*a4*n - 102.0*a4 + 24.0*a2*b2*n.pow(2) + 216.0*a2*b2*n - 528.0*a2*b2 - 16.0*b4*n.pow(2) + 48.0*b4*n - 32.0*b4)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d4(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      return a4*b1*(n + 4.0)*(3.0*a2 - b2*n + b2)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d5(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a5 = a2*a3;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      return -a5*(n + 5.0)*(a2*n - 13.0*a2 + 24.0*b2*n - 24.0*b2)/(128.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d6(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a6 = a2*a2*a2;
      const auto& b1 = this->b();
      return -a6*b1*(n + 6.0)/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y4D1::ACoeff_t I4Y4D1::d7(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a7 = a2*a2*a2*a1;
      return -a7*(n + 7.0)/(128.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   void I4Y4D1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-4, 4, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(15*std::max(this->rows(),this->cols()));
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

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
