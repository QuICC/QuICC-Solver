/**
 * @file I4Y1.cpp
 * @brief Source of the implementation of the I^4 Y sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y1.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I4Y1::I4Y1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I4Y1::ACoeff_t I4Y1::d_5(const ACoeff_t& n) const
   {
      return Internal::Math::pow(this->a()/2.0,5)/(n*(n - 3.0)*(n - 2.0)*(n - 1.0));
   }

   I4Y1::ACoeff_t I4Y1::d_4(const ACoeff_t& n) const
   {
      return Internal::Math::pow(this->a()/2.0,4)*this->b()/(n*(n - 3.0)*(n - 2.0)*(n - 1.0));
   }

   I4Y1::ACoeff_t I4Y1::d_3(const ACoeff_t& n) const
   {
      return -3.0*Internal::Math::pow(this->a()/2.0,5)/(n*(n - 2.0)*(n - 1.0)*(n + 1.0));
   }

   I4Y1::ACoeff_t I4Y1::d_2(const ACoeff_t& n) const
   {
      return -Internal::Math::pow(this->a(),4)*this->b()/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0));
   }

   I4Y1::ACoeff_t I4Y1::d_1(const ACoeff_t& n) const
   {
      return 2.0*Internal::Math::pow(this->a()/2.0,5)*(n - 8.0)/(n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0));
   }

   I4Y1::ACoeff_t I4Y1::d0(const ACoeff_t& n) const
   {
      return 6.0*Internal::Math::pow(this->a()/2.0,4)*this->b()/((n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0));
   }

   I4Y1::ACoeff_t I4Y1::d1(const ACoeff_t& n) const
   {
      return 2.0*Internal::Math::pow(this->a()/2.0,5)*(n + 8.0)/(n*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y1::ACoeff_t I4Y1::d2(const ACoeff_t& n) const
   {
      return -4.0*Internal::Math::pow(this->a()/2.0,4)*this->b()/(n*(n - 1.0)*(n + 1.0)*(n + 3.0));
   }

   I4Y1::ACoeff_t I4Y1::d3(const ACoeff_t& n) const
   {
      return -3.0*Internal::Math::pow(this->a()/2.0,5)/(n*(n - 1.0)*(n + 1.0)*(n + 2.0));
   }

   I4Y1::ACoeff_t I4Y1::d4(const ACoeff_t& n) const
   {
      return Internal::Math::pow(this->a()/2.0,4)*this->b()/(n*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y1::ACoeff_t I4Y1::d5(const ACoeff_t& n) const
   {
      return Internal::Math::pow(this->a()/2.0,5)/(n*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   void I4Y1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-4, 4, this->rows()-1);
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

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
