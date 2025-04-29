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
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

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
      using namespace Internal::Literals;
      const auto c = Internal::Math::pow(this->a()/2_mp,5);
      return c/(n*(n - 1_mp));
   }

   I2Y3::ACoeff_t I2Y3::d_4(const ACoeff_t& n) const
   {
      using namespace Internal::Literals;
      const auto& b1 = this->b();
      const auto c = Internal::Math::pow(this->a()/2_mp,4);
      return c*3_mp*b1/(n*(n - 1_mp));
   }

   I2Y3::ACoeff_t I2Y3::d_3(const ACoeff_t& n) const
   {
      using namespace Internal::Literals;
      const auto a2 = this->a()*this->a();
      const auto b2 = this->b()*this->b();
      const auto c = Internal::Math::pow(this->a()/2_mp,3);
      return c*(a2*(n + 3_mp) + 12_mp*b2*(n + 1_mp))/(4_mp*n*(n - 1_mp)*(n + 1_mp));
   }

   I2Y3::ACoeff_t I2Y3::d_2(const ACoeff_t& n) const
   {
      using namespace Internal::Literals;
      const auto& b1 = this->b();
      const auto a2 = this->a()*this->a();
      const auto b2 = this->b()*this->b();
      const auto c = Internal::Math::pow(this->a()/2_mp, 2);
      return c*b1*(3_mp*a2 + 2_mp*b2*(n + 1_mp))/(2_mp*n*(n - 1_mp)*(n + 1_mp));
   }

   I2Y3::ACoeff_t I2Y3::d_1(const ACoeff_t& n) const
   {
      using namespace Internal::Literals;
      const auto a_2 = this->a()/2_mp;
      const auto a2 = a_2*a_2;
      const auto a3 = a_2*a2;
      const auto b_2 = this->b()/2_mp;
      const auto b2 = b_2*b_2;
      return -2_mp*a3*(a2 + 6_mp*b2)/(n*(n + 1_mp));
   }

   I2Y3::ACoeff_t I2Y3::d0(const ACoeff_t& n) const
   {
      using namespace Internal::Literals;
      const auto& b1 = this->b();
      const auto a2 = this->a()*this->a();
      const auto b2 = this->b()*this->b();
      return -a2*b1*(3_mp*a2 + 4_mp*b2)/(8_mp*(n - 1_mp)*(n + 1_mp));
   }

   I2Y3::ACoeff_t I2Y3::d1(const ACoeff_t& n) const
   {
      using namespace Internal::Literals;
      return this->d_1(n - 1_mp);
   }

   I2Y3::ACoeff_t I2Y3::d2(const ACoeff_t& n) const
   {
      using namespace Internal::Literals;
      const auto& b1 = this->b();
      const auto a2 = this->a()*this->a();
      const auto b2 = this->b()*this->b();
      return -a2*b1*(3_mp*a2 - 2_mp*b2*(n - 1_mp))/(8_mp*n*(n - 1_mp)*(n + 1_mp));
   }

   I2Y3::ACoeff_t I2Y3::d3(const ACoeff_t& n) const
   {
      using namespace Internal::Literals;
      const auto a2 = this->a()*this->a();
      const auto b2 = this->b()*this->b();
      const auto c = Internal::Math::pow(this->a()/2_mp,3);
      return c*(a2*(n - 3_mp) + 12_mp*b2*(n - 1_mp))/(4_mp*n*(n - 1_mp)*(n + 1_mp));
   }

   I2Y3::ACoeff_t I2Y3::d4(const ACoeff_t& n) const
   {
      using namespace Internal::Literals;
      const auto& b1 = this->b();
      const auto c = Internal::Math::pow(this->a()/2_mp,4);
      return c*3_mp*b1/(n*(n + 1_mp));
   }

   I2Y3::ACoeff_t I2Y3::d5(const ACoeff_t& n) const
   {
      using namespace Internal::Literals;
      const auto c = Internal::Math::pow(this->a()/2_mp,5);
      return c/(n*(n + 1_mp));
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
