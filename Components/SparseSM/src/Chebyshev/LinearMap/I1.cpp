/** 
 * @file I1.cpp
 * @brief Source of the implementation of the I sparse operator, wit y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I1.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I1::I1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I1::~I1()
   {
   }

   I1::ACoeff_t I1::d_1(const ACoeff_t& n) const
   {
      return this->a()/(2.0*n);
   }

   I1::ACoeff_t I1::d1(const ACoeff_t& n) const
   {
      return -d_1(n);
   }

   void I1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-1, 1, this->rows()-1);
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
