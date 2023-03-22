/** 
 * @file ValueD2.cpp
 * @brief Source of the implementation of the boundary value and second derivative stencil
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/ValueD2.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Stencil {

   ValueD2::ValueD2(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   ValueD2::ACoeff_t ValueD2::d_4(const ACoeff_t& n) const
   {
      ACoeff_t val_num = (n - 3.0)*(2.0*n.pow(2) - 12.0*n + 19.0);
      ACoeff_t val_den = (n - 1.0)*(2.0*n.pow(2) - 4.0*n + 3.0);
      ACoeff_t val = val_num/val_den;
      val(0) = 1.0/38.0;

      return val;
   }

   ValueD2::ACoeff_t ValueD2::d_2(const ACoeff_t& n) const
   {
      ACoeff_t val_num = -2.0*n*(2.0*n.pow(2) + 7.0);
      ACoeff_t val_den = (n + 1.0)*(2.0*n.pow(2) + 4.0*n + 3.0);
      ACoeff_t val = val_num/val_den;
      val(0) = -10.0/19.0;

      return val;
   }

   ValueD2::ACoeff_t ValueD2::d0(const ACoeff_t& n) const
   {
      return ACoeff_t::Ones(n.size());
   }

   void ValueD2::buildTriplets(TripletList_t& list) const
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
