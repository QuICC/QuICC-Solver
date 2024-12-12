/**
 * @file I2Y4SphLapl.cpp
 * @brief Source of the implementation of the I^2 Y^4 spherical laplacian sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y4SphLapl.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I2Y4SphLapl::I2Y4SphLapl(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t l)
      : ISphericalOperator(rows, cols, lower, upper, l)
   {
   }

   I2Y4SphLapl::ACoeff_t I2Y4SphLapl::d_4(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      return -(a4*(l1 - n + 4.0)*(l1 + n - 3.0))/(16.0*n*(n - 1.0));
   }


   I2Y4SphLapl::ACoeff_t I2Y4SphLapl::d_3(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto& b1 = this->b();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      return -(a3*b1*(l2 + l1 - 2.0*n*n + 11.0*n - 15.0))/(4.0*n*(n - 1.0));
   }

   I2Y4SphLapl::ACoeff_t I2Y4SphLapl::d_2(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto n2 = n*n;
      const auto n3 = n2*n;
      return -a2*(a2*l2 + a2*l1 - 2*a2*n3 + 6*a2*n2 + 2*a2*n - 12*a2 + 2*b2*l2*n + 2*b2*l2 + 2*b2*l1*n + 2*b2*l1 - 12*b2*n3 + 36*b2*n2 - 48*b2)/(8*n*(n - 1)*(n + 1));
   }

   I2Y4SphLapl::ACoeff_t I2Y4SphLapl::d_1(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto n2 = n*n;
      return a1*b1*(a2*l2 + a2*l1 + 6*a2*n2 - 3*a2*n - 15*a2 + 8*b2*n2 - 4*b2*n - 12*b2)/(4*n*(n + 1));
   }

   I2Y4SphLapl::ACoeff_t I2Y4SphLapl::d0(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto b4 = b2*b2;
      const auto n2 = n*n;
      return (a4*l2 + a4*l1 + 3*a4*n2 - 9*a4 + 4*a2*b2*l2 + 4*a2*b2*l1 + 24*a2*b2*n2 - 48*a2*b2 + 8*b4*n2 - 8*b4)/(8*(n - 1)*(n + 1));
   }

   I2Y4SphLapl::ACoeff_t I2Y4SphLapl::d1(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto n2 = n*n;
      return a1*b1*(a2*l2 + a2*l1 + 6*a2*n2 + 3*a2*n - 15*a2 + 8*b2*n2 + 4*b2*n - 12*b2)/(4*n*(n - 1));
   }

   I2Y4SphLapl::ACoeff_t I2Y4SphLapl::d2(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto n2 = n*n;
      const auto n3 = n*n2;
      return a2*(a2*l2 + a2*l1 + 2*a2*n3 + 6*a2*n2 - 2*a2*n - 12*a2 - 2*b2*l2*n + 2*b2*l2 - 2*b2*l1*n + 2*b2*l1 + 12*b2*n3 + 36*b2*n2 - 48*b2)/(8*n*(n - 1)*(n + 1));
   }

   I2Y4SphLapl::ACoeff_t I2Y4SphLapl::d3(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a3 = a1*a2;
      const auto& b1 = this->b();
      const auto n2 = n*n;
      return -a3*b1*(l2 + l1 - 2*n2 - 11*n - 15)/(4*n*(n + 1));
   }

   I2Y4SphLapl::ACoeff_t I2Y4SphLapl::d4(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto a4 = a2*a2;
      return -a4*(l1 - n - 3)*(l1 + n + 4)/(16*n*(n + 1));
   }

   void I2Y4SphLapl::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-2, 2, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(13*std::max(this->rows(),this->cols()));
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
