/** 
 * @file I2Y3SphLapl.cpp
 * @brief Source of the implementation of the I^2 Y^3 spherical laplacian sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y3SphLapl.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I2Y3SphLapl::I2Y3SphLapl(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t l)
      : ISphericalOperator(rows, cols, lower, upper, l)
   {
   }

   I2Y3SphLapl::ACoeff_t I2Y3SphLapl::d_3(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto c = precision::pow(this->a()/2.0,3);
      return -c*(l1 - n + 3.0)*(l1 + n - 2.0)/(n*(n - 1.0));
   }

   I2Y3SphLapl::ACoeff_t I2Y3SphLapl::d_2(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& b1 = this->b();
      const auto c = precision::pow(this->a()/2.0,2);
      return -c*b1*(l2 + l1 - 3.0*n.pow(2) + 11.0*n - 10.0)/(n*(n - 1.0));
   }

   I2Y3SphLapl::ACoeff_t I2Y3SphLapl::d_1(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto b2 = this->b()*this->b();
      return a1*(a2*l2 + a2*l1 + 3.0*a2*n.pow(2) - a2*n - 6.0*a2 + 12.0*b2*n.pow(2) - 4.0*b2*n - 16.0*b2)/(8.0*n*(n + 1.0));
   }

   I2Y3SphLapl::ACoeff_t I2Y3SphLapl::d0(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& b1 = this->b();
      const auto a2 = this->a()*this->a();
      const auto b2 = b1*b1;
      return b1*(a2*l2 + a2*l1 + 3.0*a2*n.pow(2) - 5.0*a2 + 2.0*b2*n.pow(2) - 2.0*b2)/(2.0*(n - 1.0)*(n + 1.0));
   }

   I2Y3SphLapl::ACoeff_t I2Y3SphLapl::d1(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto b2 = this->b()*this->b();
      return a1*(a2*l2 + a2*l1 + 3.0*a2*n.pow(2) + a2*n - 6.0*a2 + 12.0*b2*n.pow(2) + 4.0*b2*n - 16.0*b2)/(8.0*n*(n - 1.0));
   }

   I2Y3SphLapl::ACoeff_t I2Y3SphLapl::d2(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto l2 = l1*l1;
      const auto& b1 = this->b();
      const auto c = precision::pow(this->a()/2.0,2);
      return -c*b1*(l2 + l1 - 3.0*n.pow(2) - 11.0*n - 10.0)/(n*(n + 1.0));
   }

   I2Y3SphLapl::ACoeff_t I2Y3SphLapl::d3(const ACoeff_t& n) const
   {
      const auto& l1 = this->l();
      const auto c = precision::pow(this->a()/2.0,3);
      return -c*(l1 - n - 2.0)*(l1 + n + 3.0)/(n*(n + 1.0));
   }

   void I2Y3SphLapl::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-2, 2, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(9*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -3, ni, this->d_3(n));
         this->convertToTriplets(list, -2, ni, this->d_2(n));
         this->convertToTriplets(list, -1, ni, this->d_1(n));
         this->convertToTriplets(list, 0, ni, this->d0(n));
         this->convertToTriplets(list, 1, ni, this->d1(n));
         this->convertToTriplets(list, 2, ni, this->d2(n));
         this->convertToTriplets(list, 3, ni, this->d3(n));
      }
   }

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
