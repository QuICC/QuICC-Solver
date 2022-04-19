/** 
 * @file I4Y3.cpp
 * @brief Source of the implementation of the I^4 Y^3 sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y3.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I4Y3::I4Y3(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I4Y3::~I4Y3()
   {
   }

   I4Y3::ACoeff_t I4Y3::d_7(const ACoeff_t& n) const
   {
      return precision::pow(this->a()/2.0,7)/(n*(n - 3.0)*(n - 2.0)*(n - 1.0));
   }

   I4Y3::ACoeff_t I4Y3::d_6(const ACoeff_t& n) const
   {
      return 3.0*precision::pow(this->a()/2.0,6)*this->b()/(n*(n - 3.0)*(n - 2.0)*(n - 1.0));
   }

   I4Y3::ACoeff_t I4Y3::d_5(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -precision::pow(this->a()/2.0,5)*(a2*(n - 11.0) - 12.0*b2*(n + 1.0))/(4.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0));
   }

   I4Y3::ACoeff_t I4Y3::d_4(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -precision::pow(this->a()/2.0,4)*this->b()*(3.0*a2*(n - 5.0) - 2.0*b2*(n + 1))/(2.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0));
   }

   I4Y3::ACoeff_t I4Y3::d_3(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -3.0*precision::pow(this->a()/2.0,5)*(a2*(n + 6.0) + 12.0*b2*(n + 2.0))/(4.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0));
   }

   I4Y3::ACoeff_t I4Y3::d_2(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -precision::pow(this->a()/2.0,4)*this->b()*(3.0*a2*(n + 17.0) + 16.0*b2*(n + 2.0))/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0));
   }

   I4Y3::ACoeff_t I4Y3::d_1(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return 3.0*precision::pow(this->a()/2.0,5)*(a2*(n.pow(2) - 5.0*n - 34.0) + 8.0*b2*(n.pow(2) - 5.0*n - 24.0))/(4.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y3::ACoeff_t I4Y3::d0(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return 3.0*precision::pow(this->a()/2.0,4)*this->b()*(a2*(n.pow(2) - 19.0) + 2.0*b2*(n.pow(2) - 9.0))/((n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y3::ACoeff_t I4Y3::d1(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return 3.0*precision::pow(this->a()/2.0,5)*(a2*(n.pow(2) + 5.0*n - 34.0) + 8.0*b2*(n.pow(2) + 5.0*n - 24.0))/(4.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y3::ACoeff_t I4Y3::d2(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -precision::pow(this->a()/2.0,4)*this->b()*(3.0*a2*(n - 17.0) + 16.0*b2*(n - 2.0))/(4.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0));
   }

   I4Y3::ACoeff_t I4Y3::d3(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -3.0*precision::pow(this->a()/2.0,5)*(a2*(n - 6.0) + 12.0*b2*(n - 2.0))/(4.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0));
   }

   I4Y3::ACoeff_t I4Y3::d4(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -precision::pow(this->a()/2.0,4)*this->b()*(3.0*a2*(n + 5.0) - 2.0*b2*(n - 1.0))/(2.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y3::ACoeff_t I4Y3::d5(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -precision::pow(this->a()/2.0,5)*(a2*(n + 11.0) - 12.0*b2*(n - 1.0))/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y3::ACoeff_t I4Y3::d6(const ACoeff_t& n) const
   {
      return 3.0*precision::pow(this->a()/2.0,6)*this->b()/(n*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y3::ACoeff_t I4Y3::d7(const ACoeff_t& n) const
   {
      return precision::pow(this->a()/2.0,7)/(n*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   void I4Y3::buildTriplets(TripletList_t& list) const
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

}
}
}
}
