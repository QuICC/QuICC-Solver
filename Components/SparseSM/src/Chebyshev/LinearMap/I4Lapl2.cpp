/** 
 * @file I4Lapl2.cpp
 * @brief Source of the implementation of the I^4 sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Lapl2.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I4Lapl2::I4Lapl2(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t k1, const Scalar_t k2)
      : IPlaneOperator(rows, cols, lower, upper, k1, k2)
   {
   }

   I4Lapl2::ACoeff_t I4Lapl2::d_4(const ACoeff_t& n) const
   {
      auto c = precision::pow(this->a()/2.0,4);
      auto laplh = -(this->k1()*this->k1() + this->k2()*this->k2());
      return laplh*laplh*c/(n*(n - 3.0)*(n - 2.0)*(n - 1.0));
   }

   I4Lapl2::ACoeff_t I4Lapl2::d_2(const ACoeff_t& n) const
   {
      // Horizontal part
      auto c = precision::pow(this->a()/2.0,4);
      auto laplh = -(this->k1()*this->k1() + this->k2()*this->k2());
      auto hz = -laplh*laplh*c*4.0/(n*(n - 3.0)*(n - 1.0)*(n + 1.0));

      // D^2 part
      auto a2 = this->a()*this->a();
      auto diff = 2.0*laplh*a2/(4.0*n*(n - 1.0));

      return hz + diff; 
   }

   I4Lapl2::ACoeff_t I4Lapl2::d0(const ACoeff_t& n) const
   {
      // Horizontal part 
      auto c = precision::pow(this->a()/2.0,4);
      auto laplh = -(this->k1()*this->k1() + this->k2()*this->k2());
      auto hz = laplh*laplh*c*6.0/((n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0));

      // D^2 part
      auto a2 = this->a()*this->a();
      auto diff = -2.0*laplh*a2/(2.0*(n - 1.0)*(n + 1.0));

      return hz + diff + 1.0;
   }

   I4Lapl2::ACoeff_t I4Lapl2::d2(const ACoeff_t& n) const
   {
      // Horizontal part 
      auto c = precision::pow(this->a()/2.0,4);
      auto laplh = -(this->k1()*this->k1() + this->k2()*this->k2());
      auto hz = -laplh*laplh*c*4.0/(n*(n - 1.0)*(n + 1.0)*(n + 3.0));

      // D^2 part
      auto a2 = this->a()*this->a();
      auto diff = 2.0*laplh*a2/(4.0*n*(n + 1.0));

      return hz + diff;
   }

   I4Lapl2::ACoeff_t I4Lapl2::d4(const ACoeff_t& n) const
   {
      auto c = precision::pow(this->a()/2.0,4);
      auto laplh = -(this->k1()*this->k1() + this->k2()*this->k2());
      return laplh*laplh*c/(n*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   void I4Lapl2::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-4, 4, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(5*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -4, ni, this->d_4(n));
         this->convertToTriplets(list, -2, ni, this->d_2(n));
         this->convertToTriplets(list, 0, ni, this->d0(n));
         this->convertToTriplets(list, 2, ni, this->d2(n));
         this->convertToTriplets(list, 4, ni, this->d4(n));
      }
   }

}
}
}
}
