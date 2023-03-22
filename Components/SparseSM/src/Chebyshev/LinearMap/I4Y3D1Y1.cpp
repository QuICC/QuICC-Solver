/** 
 * @file I4Y3D1Y1.cpp
 * @brief Source of the implementation of the I^4 Y^3 D Y sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y3D1Y1.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I4Y3D1Y1::I4Y3D1Y1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I4Y3D1Y1::~I4Y3D1Y1()
   {
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d_7(const ACoeff_t& n) const
   {
      return precision::pow(this->a()/2.0,7)*(n - 6.0)/(n*(n - 3.0)*(n - 2.0)*(n - 1.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d_6(const ACoeff_t& n) const
   {
      return precision::pow(this->a()/2.0,6)*this->b()*(4.0*n - 21.0)/(n*(n - 3.0)*(n - 2.0)*(n - 1.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d_5(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return precision::pow(this->a()/2.0,5)*(a2*(n.pow(2) + 7.0*n - 54.0) + 12.0*b2*(2.0*n.pow(2) - 7.0*n - 9.0))/(4.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d_4(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return precision::pow(this->a()/2.0,4)*this->b()*(3.0*a2*(7.0*n - 27.0) + 2.0*b2*(4.0*n.pow(2) - 11.0*n - 15.0))/(2.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d_3(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      auto a4 = a2*a2;
      auto b4 = b2*b2;
      return -precision::pow(this->a()/2.0,3)*(3.0*a4*(n.pow(2) - 4.0*n - 28.0) + 12.0*a2*b2*(2.0*n.pow(2) - 15.0*n - 38.0) - 16.0*b4*(n.pow(2) + 3.0*n + 2.0))/(16.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d_2(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -precision::pow(this->a()/2.0,4)*this->b()*(3.0*a2*(4.0*n.pow(2) - 3.0*n - 87.0) + 16.0*b2*(2.0*n.pow(2) - 5.0*n - 18.0))/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d_1(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      auto a4 = a2*a2;
      auto b4 = b2*b2;
      return -3.0*precision::pow(this->a()/2.0,3)*(a4*(n.pow(3) + 9.0*n.pow(2) - 24.0*n - 156.0) + 8.0*a2*b2*(2.0*n.pow(3) + 11.0*n.pow(2) - 33.0*n - 144.0) + 16.0*b4*(n.pow(3) + 2.0*n.pow(2) - 9.0*n - 18.0))/(16.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d0(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -3.0*precision::pow(this->a()/2.0,4)*this->b()*(a2*(7.0*n.pow(2) - 93.0) + 14.0*b2*(n.pow(2) - 9.0))/((n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d1(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      auto a4 = a2*a2;
      auto b4 = b2*b2;
      return 3.0*precision::pow(this->a()/2.0,3)*(a4*(n.pow(3) - 9.0*n.pow(2) - 24.0*n + 156.0) + 8.0*a2*b2*(2.0*n.pow(3) - 11.0*n.pow(2) - 33.0*n + 144.0) + 16.0*b4*(n - 3.0)*(n - 2.0)*(n + 3.0))/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d2(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return precision::pow(this->a()/2.0,4)*this->b()*(3.0*a2*(4.0*n.pow(2) + 3.0*n - 87.0) + 16.0*b2*(2.0*n.pow(2) + 5.0*n - 18.0))/(4.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d3(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      auto a4 = a2*a2;
      auto b4 = b2*b2;
      return precision::pow(this->a()/2.0,3)*(3.0*a4*(n.pow(2) + 4.0*n - 28.0) + 12.0*a2*b2*(2.0*n.pow(2) + 15.0*n - 38.0) - 16.0*b4*(n.pow(2) - 3.0*n + 2.0))/(16.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d4(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return precision::pow(this->a()/2.0,4)*this->b()*(3.0*a2*(7.0*n + 27.0) - 2.0*b2*(4.0*n.pow(2) + 11.0*n - 15.0))/(2.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d5(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      auto b2 = this->b()*this->b();
      return -precision::pow(this->a()/2.0,5)*(a2*(n.pow(2) - 7.0*n - 54.0) + 12.0*b2*(2.0*n.pow(2) + 7.0*n - 9.0))/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d6(const ACoeff_t& n) const
   {
      return -precision::pow(this->a()/2.0,6)*this->b()*(4.0*n + 21.0)/(n*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   I4Y3D1Y1::ACoeff_t I4Y3D1Y1::d7(const ACoeff_t& n) const
   {
      return -precision::pow(this->a()/2.0,7)*(n + 6.0)/(n*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   void I4Y3D1Y1::buildTriplets(TripletList_t& list) const
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
