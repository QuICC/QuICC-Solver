/** 
 * @file Value.cpp
 * @brief Source of the implementation of the boundary value stencil
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/Value.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Stencil {

   Value::Value(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   Value::ACoeff_t Value::d_2(const ACoeff_t& n) const
   {
      ACoeff_t val = -ACoeff_t::Ones(n.size());
      val(0) *= 0.5;
      return val;
   }

   Value::ACoeff_t Value::d0(const ACoeff_t& n) const
   {
      return ACoeff_t::Ones(n.size());
   }

   void Value::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows(), 0, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();
      ACoeffI ni2 = ni.bottomRows(this->rows()-2);
      ACoeff_t n2 = n.bottomRows(this->rows()-2);

      if(n.size() > 0)
      {
         list.reserve(2*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -2, ni2, this->d_2(n2));
         this->convertToTriplets(list, 0, ni, this->d0(n));
      }
   }

} // Stencil
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
