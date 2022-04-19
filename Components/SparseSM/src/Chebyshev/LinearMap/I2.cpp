/** 
 * @file I2.cpp
 * @brief Source of the implementation of the I^2 sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I2::I2(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I2::~I2()
   {
   }

   I2::ACoeff_t I2::d_2(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      return a2/(4.0*n*(n - 1.0));
   }

   I2::ACoeff_t I2::d0(const ACoeff_t& n) const
   {
      auto a2 = this->a()*this->a();
      return -a2/(2.0*(n - 1.0)*(n + 1.0));
   }

   I2::ACoeff_t I2::d2(const ACoeff_t& n) const
   {
      return d_2(n+1.0);
   }

   void I2::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-2, 2, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(3*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -2, ni, this->d_2(n));
         this->convertToTriplets(list, 0, ni, this->d0(n));
         this->convertToTriplets(list, 2, ni, this->d2(n));
      }
   }

}
}
}
}
