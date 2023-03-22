/** 
 * @file I4Y4.cpp
 * @brief Source of the implementation of the I^4 Y^4 sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4.hpp"


namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I4Y4::I4Y4(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I4Y4::ACoeff_t I4Y4::d_8(const ACoeff_t& n) const
   {
      const auto a8 = precision::pow(this->a(),8);
      return a8/(256.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0));
   }

   I4Y4::ACoeff_t I4Y4::d_7(const ACoeff_t& n) const
   {
      const auto a7 = precision::pow(this->a(),7);
      const auto& b1 = this->b();
      return (a7*b1)/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0));
   }

   I4Y4::ACoeff_t I4Y4::d_6(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a6 = a2*a2*a2;
      const auto b2 = b1*b1;
      return (3.0*a6*(2.0*b2*n + a2 + 2.0*b2))/(64.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 1.0));
   }

   I4Y4::ACoeff_t I4Y4::d_5(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a5 = a2*a3;
      const auto b2 = b1*b1;
      return -(a5*b1*(a2*n - 4.0*b2*n - 11.0*a2 - 4.0*b2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 1.0));
   }

   I4Y4::ACoeff_t I4Y4::d_4(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto b2 = b1*b1;
      const auto b4 = b2*b2;
      return -(a4*(a4*n.pow(2) - 19.0*a4 + 12.0*a2*b2*n.pow(2) - 36.0*a2*b2*n - 120.0*a2*b2 - 4.0*b4*n.pow(2) - 12.0*b4*n - 8.0*b4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4::ACoeff_t I4Y4::d_3(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a5 = a2*a3;
      const auto b2 = b1*b1;
      return -(3.0*a5*b1*(a2*n + 4.0*b2*n + 6.0*a2 + 8.0*b2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4::ACoeff_t I4Y4::d_2(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto b2 = b1*b1;
      const auto b4 = b2*b2;
      return -(a4*(9.0*a4*n + 33.0*a4 + 6.0*a2*b2*n.pow(2) + 120*a2*b2*n + 306.0*a2*b2 + 16.0*b4*n.pow(2) + 80.0*b4*n + 96.0*b4))/(64.0*n*(n - 1.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4::ACoeff_t I4Y4::d_1(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a5 = a2*a3;
      const auto b2 = b1*b1;
      return (a5*b1*(3.0*a2*n.pow(2) - 15.0*a2*n - 102.0*a2 + 8.0*b2*n.pow(2) - 40.0*b2*n - 192.0*b2))/(32.0*n*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4::ACoeff_t I4Y4::d0(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto b2 = b1*b1;
      const auto b4 = b2*b2;
      return (3.0*a4*(a4*n.pow(2) - 29.0*a4 + 16.0*a2*b2*n.pow(2) - 304.0*a2*b2 + 16.0*b4*n.pow(2) - 144.0*b4))/(128.0*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4::ACoeff_t I4Y4::d1(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a5 = a2*a3;
      const auto b2 = b1*b1;
      return (a5*b1*(3.0*a2*n.pow(2) + 15.0*a2*n - 102.0*a2 + 8.0*b2*n.pow(2) + 40.0*b2*n - 192.0*b2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0));
   }

   I4Y4::ACoeff_t I4Y4::d2(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto b2 = b1*b1;
      const auto b4 = b2*b2;
      return (a4*(9.0*a4*n - 33.0*a4 - 6.0*a2*b2*n.pow(2) + 120.0*a2*b2*n - 306.0*a2*b2 - 16.0*b4*n.pow(2) + 80.0*b4*n - 96.0*b4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 1.0));
   }

   I4Y4::ACoeff_t I4Y4::d3(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a5 = a2*a3;
      const auto b2 = b1*b1;
      return -(3.0*a5*b1*(a2*n + 4.0*b2*n - 6.0*a2 - 8.0*b2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4::ACoeff_t I4Y4::d4(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto b2 = b1*b1;
      const auto b4 = b2*b2;
      return -(a4*(a4*n.pow(2) - 19.0*a4 + 12.0*a2*b2*n.pow(2) + 36.0*a2*b2*n - 120.0*a2*b2 - 4.0*b4*n.pow(2) + 12.0*b4*n - 8.0*b4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n + 3.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4::ACoeff_t I4Y4::d5(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a5 = a2*a3;
      const auto b2 = b1*b1;
      return -(a5*b1*(a2*n - 4.0*b2*n + 11.0*a2 + 4.0*b2))/(32.0*n*(n - 1.0)*(n + 3.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4::ACoeff_t I4Y4::d6(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto a6 = a2*a4;
      const auto b2 = b1*b1;
      return -(3.0*a6*(a2 - 2.0*b2*n + 2.0*b2))/(64.0*n*(n - 1.0)*(n + 3.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4::ACoeff_t I4Y4::d7(const ACoeff_t& n) const
   {
      return this->d_7(n + 3.0);
   }

   I4Y4::ACoeff_t I4Y4::d8(const ACoeff_t& n) const
   {
      return this->d_8(n + 3.0);
   }

   void I4Y4::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-4, 4, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(15*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -8, ni, this->d_8(n));
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
         this->convertToTriplets(list, 8, ni, this->d8(n));
      }
   }

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
