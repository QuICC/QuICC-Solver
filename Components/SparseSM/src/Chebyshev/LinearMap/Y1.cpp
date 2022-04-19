/** 
 * @file Y1.cpp
 * @brief Source of the implementation of the Y sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y1.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   Y1::Y1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   Y1::~Y1()
   {
   }

   Y1::ACoeff_t Y1::d_1(const ACoeff_t& n) const
   {
      return ACoeff_t::Constant(n.size(),this->a()/2.0);
   }

   Y1::ACoeff_t Y1::d0(const ACoeff_t& n) const
   {
      return ACoeff_t::Constant(n.size(),this->b());
   }

   Y1::ACoeff_t Y1::d1(const ACoeff_t& n) const
   {
      return d_1(n);
   }

   void Y1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows(), 0, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(3*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -1, ni, this->d_1(n));
         this->convertToTriplets(list, 0, ni, this->d0(n));
         this->convertToTriplets(list, 1, ni, this->d1(n));
      }
   }

}
}
}
}
