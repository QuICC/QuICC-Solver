/** 
 * @file R1D1DivR1.cpp
 * @brief Source of the implementation of the boundary value of r D 1/r stencil
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/R1D1DivR1.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Stencil {

   R1D1DivR1::R1D1DivR1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   R1D1DivR1::ACoeff_t R1D1DivR1::d_2(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto val_num = (n - 2.0)*(b2*(n - 2.0)*(n - 1.0) - a2*(n - 3.0)*n);
      const auto val_den = n*(a2*(n - 2.0)*(n + 1.0) - b2*(n - 1.0)*n);
      ACoeff_t val = val_num/val_den;
      val(0) = 0;

      return val;
   }

   R1D1DivR1::ACoeff_t R1D1DivR1::d_1(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto val_num = -4.0*a1*b1;
      const auto val_den = (n + 1.0)*(a2*(n.pow(2) + n - 2.0) - b2*n*(n + 1.0));
      ACoeff_t val = val_num/val_den;
      val(0) = a1/(2.0*b1);

      return val;
   }

   R1D1DivR1::ACoeff_t R1D1DivR1::d0(const ACoeff_t& n) const
   {
      return ACoeff_t::Ones(n.size());
   }

   void R1D1DivR1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows(), 0, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();
      ACoeffI ni1 = ni.bottomRows(this->rows()-1);
      ACoeff_t n1 = n.bottomRows(this->rows()-1);
      ACoeffI ni2 = ni.bottomRows(this->rows()-2);
      ACoeff_t n2 = n.bottomRows(this->rows()-2);

      if(n.size() > 0)
      {
         list.reserve(2*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -2, ni2, this->d_2(n2));
         this->convertToTriplets(list, -1, ni1, this->d_1(n1));
         this->convertToTriplets(list, 0, ni, this->d0(n));
      }
   }

} // Stencil
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
