/**
 * @file I2Y4.cpp
 * @brief Source of the implementation of the I^2 Y^4 sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y4.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I2Y4::I2Y4(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I2Y4::ACoeff_t I2Y4::d_6(const ACoeff_t& n) const
   {
      auto a1 = this->a();
      return Internal::Math::pow(a1/2.0,6)/(n*(n - 1.0));
   }

   I2Y4::ACoeff_t I2Y4::d_5(const ACoeff_t& n) const
   {
      auto a1 = this->a();
      return Internal::Math::pow(a1,5)*this->b()/(8.0*n*(n - 1.0));
   }

   I2Y4::ACoeff_t I2Y4::d_4(const ACoeff_t& n) const
   {
      auto a1 = this->a();
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();

      return Internal::Math::pow(a1/2.0,4)*(a2*(n + 2.0) + 12.0*b2*(n + 1.0))/(2.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y4::ACoeff_t I2Y4::d_3(const ACoeff_t& n) const
   {
      auto a1 = this->a();
      auto b1 = this->b();
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return Internal::Math::pow(a1/2.0,3)*b1*(a2*(n + 3.0) + 4.0*b2*(n + 1.0))/(n*(n - 1.0)*(n + 1.0));
   }

   I2Y4::ACoeff_t I2Y4::d_2(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      auto a4 = a2*a2;
      auto b4 = b2*b2;

      return -a2*(a4*(n - 5.0) - 48.0*a2*b2 - 16.0*b4*(n + 1.0))/(64.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y4::ACoeff_t I2Y4::d_1(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      auto a3 = a2*this->a();
      return -a3*this->b()*(a2 + 2.0*b2)/(4.0*n*(n + 1.0));
   }

   I2Y4::ACoeff_t I2Y4::d0(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      auto a4 = a2*a2;
      auto b4 = b2*b2;
      return -a2*(a4 + 12.0*a2*b2 + 8.0*b4)/(16.0*(n - 1.0)*(n + 1.0));
   }

   I2Y4::ACoeff_t I2Y4::d1(const ACoeff_t& n) const
   {
      return this->d_1(n - 1.0);
   }

   I2Y4::ACoeff_t I2Y4::d2(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      auto a4 = a2*a2;
      auto b4 = b2*b2;
      return -a2*(a4*(n + 5.0) + 48.0*a2*b2 - 16.0*b4*(n - 1.0))/(64.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y4::ACoeff_t I2Y4::d3(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      auto b1 = this->b();
      auto a3 = a2*this->a();
      return a3*b1*(a2*(n - 3.0) + 4.0*b2*(n - 1.0))/(8.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y4::ACoeff_t I2Y4::d4(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();

      return Internal::Math::pow(this->a()/2.0,4)*(a2*(n - 2.0) + 12.0*b2*(n - 1.0))/(2.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y4::ACoeff_t I2Y4::d5(const ACoeff_t& n) const
   {
      return d_5(n + 1.0);
   }

   I2Y4::ACoeff_t I2Y4::d6(const ACoeff_t& n) const
   {
      return d_6(n + 1.0);
   }

   void I2Y4::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-2, 2, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(13*std::max(this->rows(),this->cols()));
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

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
