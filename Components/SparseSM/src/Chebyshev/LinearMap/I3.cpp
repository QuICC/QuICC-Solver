/**
 * @file I3.cpp
 * @brief Source of the implementation of the I^3 sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I3.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I3::I3(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I3::ACoeff_t I3::d_3(const ACoeff_t& n) const
   {
      return Internal::Math::pow(this->a()/2.0,3)/(n*(n - 2.0)*(n - 1.0));
   }

   I3::ACoeff_t I3::d_1(const ACoeff_t& n) const
   {
      return -Internal::Math::pow(this->a()/2.0,3)*3.0/(n*(n - 2.0)*(n + 1.0));
   }

   I3::ACoeff_t I3::d1(const ACoeff_t& n) const
   {
      return Internal::Math::pow(this->a()/2.0,3)*3.0/(n*(n - 1.0)*(n + 2.0));
   }

   I3::ACoeff_t I3::d3(const ACoeff_t& n) const
   {
      return -Internal::Math::pow(this->a()/2.0,3)/(n*(n + 1.0)*(n + 2.0));
   }

   void I3::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-3, 3, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(4*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -3, ni, this->d_3(n));
         this->convertToTriplets(list, -1, ni, this->d_1(n));
         this->convertToTriplets(list, 1, ni, this->d1(n));
         this->convertToTriplets(list, 3, ni, this->d3(n));
      }
   }

}
}
}
}
