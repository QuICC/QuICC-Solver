/** 
 * @file D1.cpp
 * @brief Source of the implementation of the boundary value stencil
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/D1.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Stencil {

   D1::D1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   D1::ACoeff_t D1::d_2(const ACoeff_t& n) const
   {
      ACoeff_t val = -(n - 2.0).pow(2)/n.pow(2);
      return val;
   }

   D1::ACoeff_t D1::d0(const ACoeff_t& n) const
   {
      return ACoeff_t::Ones(n.size());
   }

   void D1::buildTriplets(TripletList_t& list) const
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
