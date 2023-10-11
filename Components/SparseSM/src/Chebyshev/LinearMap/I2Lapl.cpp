/**
 * @file I2Lapl.cpp
 * @brief Source of the implementation of the I^2 plane layer Cartesian laplacian sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Lapl.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   I2Lapl::I2Lapl(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t k1, const Scalar_t k2)
      : IPlaneOperator(rows, cols, lower, upper, k1, k2)
   {
   }

   I2Lapl::ACoeff_t I2Lapl::d_2(const ACoeff_t& n) const
   {
      auto laplh = -(this->k1()*this->k1() + this->k2()*this->k2());
      auto a2 = this->a()*this->a();
      return a2*laplh/(4.0*n*(n - 1.0));
   }

   I2Lapl::ACoeff_t I2Lapl::d0(const ACoeff_t& n) const
   {
      auto laplh = -(this->k1()*this->k1() + this->k2()*this->k2());
      auto a2 = this->a()*this->a();
      return -a2*laplh/(2.0*(n - 1.0)*(n + 1.0)) + 1.0;
   }

   I2Lapl::ACoeff_t I2Lapl::d2(const ACoeff_t& n) const
   {
      return d_2(n+1.0);
   }

   void I2Lapl::buildTriplets(TripletList_t& list) const
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
