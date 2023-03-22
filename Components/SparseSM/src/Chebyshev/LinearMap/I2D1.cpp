/** 
 * @file I2D1.cpp
 * @brief Source of the implementation of the I^2 D sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2D1.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I2D1::I2D1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I2D1::~I2D1()
   {
   }

   I2D1::ACoeff_t I2D1::d_1(const ACoeff_t& n) const
   {
      return this->a()/(2.0*n);
   }

   I2D1::ACoeff_t I2D1::d1(const ACoeff_t& n) const
   {
      return -d_1(n);
   }

   void I2D1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-2, 2, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(2*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -1, ni, this->d_1(n));
         this->convertToTriplets(list, 1, ni, this->d1(n));
      }
   }

}
}
}
}
