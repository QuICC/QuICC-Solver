/** 
 * @file I4.cpp
 * @brief Source of the implementation of the I^4 sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I4::I4(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I4::ACoeff_t I4::d_4(const ACoeff_t& n) const
   {
      return precision::pow(this->a()/2.0,4)/(n*(n - 3.0)*(n - 2.0)*(n - 1.0));
   }

   I4::ACoeff_t I4::d_2(const ACoeff_t& n) const
   {
      return -precision::pow(this->a()/2.0,4)*4.0/(n*(n - 3.0)*(n - 1.0)*(n + 1.0)); 
   }

   I4::ACoeff_t I4::d0(const ACoeff_t& n) const
   {
      return precision::pow(this->a()/2.0,4)*6.0/((n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0));
   }

   I4::ACoeff_t I4::d2(const ACoeff_t& n) const
   {
      return -precision::pow(this->a()/2.0,4)*4.0/(n*(n - 1.0)*(n + 1.0)*(n + 3.0));
   }

   I4::ACoeff_t I4::d4(const ACoeff_t& n) const
   {
      return precision::pow(this->a()/2.0,4)/(n*(n + 1.0)*(n + 2.0)*(n + 3.0));
   }

   void I4::buildTriplets(TripletList_t& list) const
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
