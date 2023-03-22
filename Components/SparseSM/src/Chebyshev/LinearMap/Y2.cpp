/** 
 * @file Y2.cpp
 * @brief Source of the implementation of the Y^2 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y2.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   Y2::Y2(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   Y2::~Y2()
   {
   }

   Y2::ACoeff_t Y2::d_2(const ACoeff_t& n) const
   {
      return ACoeff_t::Constant(n.size(), precision::pow(this->a(),2)/4.0);
   }

   Y2::ACoeff_t Y2::d_1(const ACoeff_t& n) const
   {
      return ACoeff_t::Constant(n.size(), this->a()*this->b());
   }

   Y2::ACoeff_t Y2::d0(const ACoeff_t& n) const
   {
      return ACoeff_t::Constant(n.size(), (precision::pow(this->a(),2) + 2.0*precision::pow(this->b(),2))/2.0);
   }

   Y2::ACoeff_t Y2::d1(const ACoeff_t& n) const
   {
      return d_1(n);
   }

   Y2::ACoeff_t Y2::d2(const ACoeff_t& n) const
   {
      return d_2(n);
   }

   void Y2::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows(), 0, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(5*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -2, ni, this->d_2(n));
         this->convertToTriplets(list, -1, ni, this->d_1(n));
         this->convertToTriplets(list, 0, ni, this->d0(n));
         this->convertToTriplets(list, 1, ni, this->d1(n));
         this->convertToTriplets(list, 2, ni, this->d2(n));
      }
   }

}
}
}
}
