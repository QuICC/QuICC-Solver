/** 
 * @file I4Y4SphLapl.cpp
 * @brief Source of the implementation of the I^4 Y^4 spherical laplacian sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4SphLapl.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I4Y4SphLapl::I4Y4SphLapl(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t l)
      : ISphericalOperator(rows, cols, lower, upper, l)
   {
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d_6(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto a6 = a2*a4;
      return -(a6*(l1 - n + 6.0)*(l1 + n - 5.0))/(64.0*n*(n - 3.0)*(n - 1.0)*(n - 2.0));
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d_5(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a5 = a2*a3;
      return -(a5*b1*(l2 + l1 - 2.0*n.pow(2) + 19.0*n - 45.0))/(16.0*n*(n - 3.0)*(n - 1.0)*(n - 2.0));
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d_4(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      const auto a4 = a2*a2;
      return (a4*(a2*l2*n - 5.0*a2*l2 + a2*l1*n - 5.0*a2*l1 + a2*n.pow(3) - 3.0*a2*n.pow(2) - 28.0*a2*n + 96.0*a2 - 2.0*b2*l2*n - 2.0*b2*l2 - 2.0*b2*l1*n - 2.0*b2*l1 + 12.0*b2*n.pow(3) - 84.0*b2*n.pow(2) + 96.0*b2*n + 192.0*b2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 1.0));
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d_3(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      const auto a3 = a1*a2;
      return (a3*b1*(3.0*a2*l2 + 3.0*a2*l1 + 2.0*a2*n.pow(2) + 11.0*a2*n - 75.0*a2 + 8.0*b2*n.pow(2) - 20.0*b2*n - 28.0*b2))/(16.0*n*(n - 1.0)*(n - 2.0)*(n + 1.0));
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d_2(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      const auto a4 = a2*a2;
      const auto b4 = b2*b2;
      return (a2*(a4*l2*n + 17.0*a4*l2 + a4*l1*n + 17.0*a4*l1 - a4*n.pow(3) + 24.0*a4*n.pow(2) - 5.0*a4*n - 294.0*a4 + 16.0*a2*b2*l2*n + 32.0*a2*b2*l2 + 16.0*a2*b2*l1*n + 32.0*a2*b2*l1 + 192.0*a2*b2*n.pow(2) - 288.0*a2*b2*n - 1344.0*a2*b2 + 16.0*b4*n.pow(3) - 112.0*b4*n - 96.0*b4))/(64.0*n*(n - 1.0)*(n - 3.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d_1(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      const auto a3 = a1*a2;
      return -(a3*b1*(a2*l2*n - 8.0*a2*l2 + a2*l1*n - 8.0*a2*l1 + 2.0*a2*n.pow(3) - 15.0*a2*n.pow(2) - 23.0*a2*n + 180.0*a2 + 4.0*b2*n.pow(3) - 30.0*b2*n.pow(2) + 2.0*b2*n + 156.0*b2))/(8.0*n*(n - 2.0)*(n - 3.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d0(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      const auto a4 = a2*a2;
      const auto b4 = b2*b2;
      return -(a2*(a4*l2*n.pow(2) - 19.0*a4*l2 + a4*l1*n.pow(2) - 19.0*a4*l1 + a4*n.pow(4) - 37.0*a4*n.pow(2) + 312.0*a4 + 6.0*a2*b2*l2*n.pow(2) - 54.0*a2*b2*l2 + 6.0*a2*b2*l1*n.pow(2) - 54.0*a2*b2*l1 + 12.0*a2*b2*n.pow(4) - 300.0*a2*b2*n.pow(2) + 1728.0*a2*b2 + 8.0*b4*n.pow(4) - 104.0*b4*n.pow(2) + 288.0*b4))/(16.0*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d1(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      const auto a3 = a1*a2;
      return -(a3*b1*(a2*l2*n + 8.0*a2*l2 + a2*l1*n + 8.0*a2*l1 + 2.0*a2*n.pow(3) + 15.0*a2*n.pow(2) - 23.0*a2*n - 180.0*a2 + 4.0*b2*n.pow(3) + 30.0*b2*n.pow(2) + 2.0*b2*n - 156.0*b2))/(8.0*n*(n - 1.0)*(n - 2.0)*(n + 3.0)*(n + 2.0));
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d2(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      const auto a4 = a2*a2;
      const auto b4 = b2*b2;
      return (a2*(a4*l2*n - 17.0*a4*l2 + a4*l1*n - 17.0*a4*l1 - a4*n.pow(3) - 24.0*a4*n.pow(2) - 5.0*a4*n + 294.0*a4 + 16.0*a2*b2*l2*n - 32.0*a2*b2*l2 + 16.0*a2*b2*l1*n - 32.0*a2*b2*l1 - 192.0*a2*b2*n.pow(2) - 288.0*a2*b2*n + 1344.0*a2*b2 + 16.0*b4*n.pow(3) - 112.0*b4*n + 96.0*b4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n + 3.0)*(n + 1.0));
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d3(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      const auto a3 = a1*a2;
      return (a3*b1*(3.0*a2*l2 + 3.0*a2*l1 + 2.0*a2*n.pow(2) - 11.0*a2*n - 75.0*a2 + 8.0*b2*n.pow(2) + 20.0*b2*n - 28.0*b2))/(16.0*n*(n - 1.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d4(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      const auto a4 = a2*a2;
      return (a4*(a2*l2*n + 5.0*a2*l2 + a2*l1*n + 5.0*a2*l1 + a2*n.pow(3) + 3.0*a2*n.pow(2) - 28.0*a2*n - 96.0*a2 - 2.0*b2*l2*n + 2.0*b2*l2 - 2.0*b2*l1*n + 2.0*b2*l1 + 12.0*b2*n.pow(3) + 84.0*b2*n.pow(2) + 96.0*b2*n - 192.0*b2))/(32.0*n*(n - 1.0)*(n + 3.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d5(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a5 = a2*a3;
      return -(a5*b1*(l2 + l1 - 2.0*n.pow(2) - 19.0*n - 45.0))/(16.0*n*(n + 3.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4SphLapl::ACoeff_t I4Y4SphLapl::d6(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto a6 = a3*a3;
      return -(a6*(l1 - n - 5.0)*(l1 + n + 6.0))/(64.0*n*(n + 3.0)*(n + 2.0)*(n + 1.0));
   }

   void I4Y4SphLapl::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-4, 4, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(9*std::max(this->rows(),this->cols()));
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
