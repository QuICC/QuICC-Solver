/**
 * @file I4Y4SphLapl2.cpp
 * @brief Source of the implementation of the I^4 Y^4 spherical bilaplacian sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4SphLapl2.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I4Y4SphLapl2::I4Y4SphLapl2(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t l)
      : ISphericalOperator(rows, cols, lower, upper, l)
   {
   }

   I4Y4SphLapl2::ACoeff_t I4Y4SphLapl2::d_4(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto c = Internal::Math::pow(this->a()/2.0,4);
      return (c*(l1 - n + 6.0)*(l1 + n - 5.0)*(l1 - n + 4.0)*(l1 + n - 3.0))/(n*(n - 3.0)*(n - 1.0)*(n - 2.0));
   }

   I4Y4SphLapl2::ACoeff_t I4Y4SphLapl2::d_3(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& b1 = this->b();
      const auto a3 = this->a()*this->a()*this->a();
      return -(a3*b1*(n - 4.0)*(l2 + l1 - n.pow(2) + 8.0*n - 15.0))/(2.0*n*(n - 1.0)*(n - 2.0));
   }

   I4Y4SphLapl2::ACoeff_t I4Y4SphLapl2::d_2(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto l3 = l1*l2;
      const auto l4 = l2*l2;
      const auto a2 = this->a()*this->a();
      const auto b2 = this->b()*this->b();
      return -(a2*(a2*l4 + 2.0*a2*l3 + 5.0*a2*l2*n - 20.0*a2*l2 + 5.0*a2*l1*n - 21.0*a2*l1 - a2*n.pow(4) + 9.0*a2*n.pow(3) - 17.0*a2*n.pow(2) - 39.0*a2*n + 108.0*a2 + 2.0*b2*l2*n.pow(2) - 4.0*b2*l2*n - 6.0*b2*l2 + 2.0*b2*l1*n.pow(2) - 4.0*b2*l1*n - 6.0*b2*l1 - 6.0*b2*n.pow(4) + 54.0*b2*n.pow(3) - 138.0*b2*n.pow(2) + 18.0*b2*n + 216.0*b2))/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0));
   }

   I4Y4SphLapl2::ACoeff_t I4Y4SphLapl2::d_1(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      return (a1*b1*(a2*l2*n - 8.0*a2*l2 + a2*l1*n - 8.0*a2*l1 + 3.0*a2*n.pow(3) - 12.0*a2*n.pow(2) - 15.0*a2*n + 72.0*a2 + 4.0*b2*n.pow(3) - 16.0*b2*n.pow(2) + 4.0*b2*n + 24.0*b2))/(2.0*n*(n + 1.0)*(n - 2.0));
   }

   I4Y4SphLapl2::ACoeff_t I4Y4SphLapl2::d0(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto l3 = l1*l2;
      const auto l4 = l2*l2;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      const auto a4 = a2*a2;
      const auto b4 = b2*b2;
      return (3.0*a4*l4 + 6.0*a4*l3 + 2.0*a4*l2*n.pow(2) - 47.0*a4*l2 + 2.0*a4*l1*n.pow(2) - 50.0*a4*l1 + 3.0*a4*n.pow(4) - 51.0*a4*n.pow(2) + 228.0*a4 + 8.0*a2*b2*l2*n.pow(2) - 32.0*a2*b2*l2 + 8.0*a2*b2*l1*n.pow(2) - 32.0*a2*b2*l1 + 24.0*a2*b2*n.pow(4) - 264.0*a2*b2*n.pow(2) + 672.0*a2*b2 + 8.0*b4*n.pow(4) - 40.0*b4*n.pow(2) + 32.0*b4)/(8.0*(n - 1.0)*(n - 2.0)*(n + 2.0)*(n + 1.0));
   }

   I4Y4SphLapl2::ACoeff_t I4Y4SphLapl2::d1(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      return (a1*b1*(a2*l2*n + 8.0*a2*l2 + a2*l1*n + 8.0*a2*l1 + 3.0*a2*n.pow(3) + 12.0*a2*n.pow(2) - 15.0*a2*n - 72.0*a2 + 4.0*b2*n.pow(3) + 16.0*b2*n.pow(2) + 4.0*b2*n - 24.0*b2))/(2.0*n*(n + 2.0)*(n - 1.0));
   }

   I4Y4SphLapl2::ACoeff_t I4Y4SphLapl2::d2(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto l3 = l1*l2;
      const auto l4 = l2*l2;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto b2 = b1*b1;
      return -(a2*(a2*l4 + 2.0*a2*l3 - 5.0*a2*l2*n - 20.0*a2*l2 - 5.0*a2*l1*n - 21.0*a2*l1 - a2*n.pow(4) - 9.0*a2*n.pow(3) - 17.0*a2*n.pow(2) + 39.0*a2*n + 108.0*a2 + 2.0*b2*l2*n.pow(2) + 4.0*b2*l2*n - 6.0*b2*l2 + 2.0*b2*l1*n.pow(2) + 4.0*b2*l1*n - 6.0*b2*l1 - 6.0*b2*n.pow(4) - 54.0*b2*n.pow(3) - 138.0*b2*n.pow(2) - 18.0*b2*n + 216.0*b2))/(4.0*n*(n + 3.0)*(n - 1.0)*(n + 1.0));
   }

   I4Y4SphLapl2::ACoeff_t I4Y4SphLapl2::d3(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      return -(a3*b1*(n + 4.0)*(l2 + l1 - n.pow(2) - 8.0*n - 15.0))/(2.0*n*(n + 2.0)*(n + 1.0));
   }

   I4Y4SphLapl2::ACoeff_t I4Y4SphLapl2::d4(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      return (a4*(l1 + n + 6.0)*(l1 - n - 3.0)*(l1 + n + 4.0)*(l1 - n - 5.0))/(16.0*n*(n + 3.0)*(n + 2.0)*(n + 1.0));
   }

   void I4Y4SphLapl2::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-4, 4, this->rows()-1);
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
