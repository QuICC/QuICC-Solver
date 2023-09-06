/** 
 * @file I2Y2.cpp
 * @brief Source of the implementation of the I^2 Y^2 sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I2Y2::I2Y2(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I2Y2::ACoeff_t I2Y2::d_4(const ACoeff_t& n) const
   {
      return precision::pow(this->a()/2.0,4)/(n*(n - 1.0));
   }

   I2Y2::ACoeff_t I2Y2::d_3(const ACoeff_t& n) const
   {
      return precision::pow(this->a(),3)*this->b()/(4.0*n*(n - 1.0));
   }

   I2Y2::ACoeff_t I2Y2::d_2(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return a2*(a2 + 2.0*b2*(n + 1.0))/(8.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y2::ACoeff_t I2Y2::d_1(const ACoeff_t& n) const
   {
      return -this->d_3(n + 1.0);
   }

   I2Y2::ACoeff_t I2Y2::d0(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -a2*(a2 + 4.0*b2)/(8.0*(n - 1.0)*(n + 1.0));
   }

   I2Y2::ACoeff_t I2Y2::d1(const ACoeff_t& n) const
   {
      return this->d_1(n - 1.0);
   }

   I2Y2::ACoeff_t I2Y2::d2(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -a2*(a2 - 2.0*b2*(n - 1.0))/(8.0*n*(n - 1.0)*(n + 1.0));
   }

   I2Y2::ACoeff_t I2Y2::d3(const ACoeff_t& n) const
   {
      return d_3(n + 1.0);
   }

   I2Y2::ACoeff_t I2Y2::d4(const ACoeff_t& n) const
   {
      return d_4(n + 1.0);
   }

   void I2Y2::buildTriplets(TripletList_t& list) const
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

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
