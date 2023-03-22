/** 
 * @file ValueD1.cpp
 * @brief Source of the implementation of the boundary value and first derivative stencil
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/ValueD1.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Stencil {

   ValueD1::ValueD1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   ValueD1::ACoeff_t ValueD1::d_4(const ACoeff_t& n) const
   {
      ACoeff_t val = (n - 3.0)/(n - 1.0);
      val(0) = 1.0/6.0;
      return val;
   }

   ValueD1::ACoeff_t ValueD1::d_2(const ACoeff_t& n) const
   {
      ACoeff_t val = -2.0*n/(n + 1.0);
      val(0) = -2.0/3.0;
      return val;
   }

   ValueD1::ACoeff_t ValueD1::d0(const ACoeff_t& n) const
   {
      return ACoeff_t::Ones(n.size());
   }

   void ValueD1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows(), 0, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();
      ACoeffI ni2 = ni.bottomRows(this->rows()-2);
      ACoeff_t n2 = n.bottomRows(this->rows()-2);
      ACoeffI ni4 = ni.bottomRows(this->rows()-4);
      ACoeff_t n4 = n.bottomRows(this->rows()-4);

      if(n.size() > 0)
      {
         list.reserve(3*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -4, ni4, this->d_4(n4));
         this->convertToTriplets(list, -2, ni2, this->d_2(n2));
         this->convertToTriplets(list, 0, ni, this->d0(n));
      }
   }

} // Stencil
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
