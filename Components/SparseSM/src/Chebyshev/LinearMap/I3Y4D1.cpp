/**
 * @file I3Y4D1.cpp
 * @brief Source of the implementation of the I^3 Y^4 D sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I3Y4D1.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   using namespace Internal::Literals;

   I3Y4D1::I3Y4D1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : ILinearMapOperator(rows, cols, lower, upper)
   {
   }

   I3Y4D1::ACoeff_t I3Y4D1::d_6(const ACoeff_t& n) const
   {
      return a<6>()*(n - 6.0_mp)/(64.0_mp*n*(n - 2.0_mp)*(n - 1.0_mp));
   }

   I3Y4D1::ACoeff_t I3Y4D1::d_5(const ACoeff_t& n) const
   {
      return a<5>()*b<1>()*(n - 5.0_mp)/(8.0_mp*n*(n - 2.0_mp)*(n - 1.0_mp));
   }

   I3Y4D1::ACoeff_t I3Y4D1::d_4(const ACoeff_t& n) const
   {
      return a<4>()*(n - 4.0_mp)*(a<2>()*n + 4.0_mp*a<2>() + 12.0_mp*b<2>()*n + 12.0_mp*b<2>())/(32.0_mp*n*(n - 2.0_mp)*(n - 1.0_mp)*(n + 1.0_mp));
   }

   I3Y4D1::ACoeff_t I3Y4D1::d_3(const ACoeff_t& n) const
   {
      return a<3>()*b<1>()*(n - 3.0_mp)*(a<2>()*n + 7.0_mp*a<2>() + 4.0_mp*b<2>()*n + 4.0_mp*b<2>())/(8.0_mp*n*(n - 2.0_mp)*(n - 1.0_mp)*(n + 1.0_mp));
   }

   I3Y4D1::ACoeff_t I3Y4D1::d_2(const ACoeff_t& n) const
   {
      return -a<2>()*(a<4>()*n*n - 15.0_mp*a<4>()*n - 46.0_mp*a<4>() - 144.0_mp*a<2>()*b<2>()*n - 288.0_mp*a<2>()*b<2>() - 16.0_mp*b<4>()*n*n - 48.0_mp*b<4>()*n - 32.0_mp*b<4>())/(64.0_mp*n*(n - 1.0_mp)*(n + 1.0_mp)*(n + 2.0_mp));
   }

   I3Y4D1::ACoeff_t I3Y4D1::d_1(const ACoeff_t& n) const
   {
      return -a<3>()*b<1>()*(a<2>()*n*n - 3.0_mp*a<2>()*n - 16.0_mp*a<2>() + 2.0_mp*b<2>()*n*n - 6.0_mp*b<2>()*n - 20.0_mp*b<2>())/(4.0_mp*n*(n - 2.0_mp)*(n + 1.0_mp)*(n + 2.0_mp));
   }

   I3Y4D1::ACoeff_t I3Y4D1::d0(const ACoeff_t& n) const
   {
      return -a<2>()*(a<4>()*n*n - 16.0_mp*a<4>() + 12.0_mp*a<2>()*b<2>()*n*n - 120.0_mp*a<2>()*b<2>() + 8.0_mp*b<4>()*n*n - 32.0_mp*b<4>())/(16.0_mp*(n - 2.0_mp)*(n - 1.0_mp)*(n + 1.0_mp)*(n + 2.0_mp));
   }

   I3Y4D1::ACoeff_t I3Y4D1::d1(const ACoeff_t& n) const
   {
      return -a<3>()*b<1>()*(a<2>()*n*n + 3.0_mp*a<2>()*n - 16.0_mp*a<2>() + 2.0_mp*b<2>()*n*n + 6.0_mp*b<2>()*n - 20.0_mp*b<2>())/(4.0_mp*n*(n - 2.0_mp)*(n - 1.0_mp)*(n + 2.0_mp));
   }

   I3Y4D1::ACoeff_t I3Y4D1::d2(const ACoeff_t& n) const
   {
      return -a<2>()*(a<4>()*n*n + 15.0_mp*a<4>()*n - 46.0_mp*a<4>() + 144.0_mp*a<2>()*b<2>()*n - 288.0_mp*a<2>()*b<2>() - 16.0_mp*b<4>()*n*n + 48.0_mp*b<4>()*n - 32.0_mp*b<4>())/(64.0_mp*n*(n - 2.0_mp)*(n - 1.0_mp)*(n + 1.0_mp));
   }

   I3Y4D1::ACoeff_t I3Y4D1::d3(const ACoeff_t& n) const
   {
      return a<3>()*b<1>()*(n + 3.0_mp)*(a<2>()*n - 7.0_mp*a<2>() + 4.0_mp*b<2>()*n - 4.0_mp*b<2>())/(8.0_mp*n*(n - 1.0_mp)*(n + 1.0_mp)*(n + 2.0_mp));
   }

   I3Y4D1::ACoeff_t I3Y4D1::d4(const ACoeff_t& n) const
   {
      return a<4>()*(n + 4.0_mp)*(a<2>()*n - 4.0_mp*a<2>() + 12.0_mp*b<2>()*n - 12.0_mp*b<2>())/(32.0_mp*n*(n - 1.0_mp)*(n + 1.0_mp)*(n + 2.0_mp));
   }

   I3Y4D1::ACoeff_t I3Y4D1::d5(const ACoeff_t& n) const
   {
      return a<5>()*b<1>()*(n + 5.0_mp)/(8.0_mp*n*(n + 1.0_mp)*(n + 2.0_mp));
   }

   I3Y4D1::ACoeff_t I3Y4D1::d6(const ACoeff_t& n) const
   {
      return a<6>()*(n + 6.0_mp)/(64.0_mp*n*(n + 1.0_mp)*(n + 2.0_mp));
   }

   void I3Y4D1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-3, 3, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(13*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -6, ni, this->d_6(n));
         this->convertToTriplets(list, -5, ni, this->d_5(n));
         this->convertToTriplets(list, -4, ni, this->d_4(n));
         this->convertToTriplets(list, -3, ni, this->d_3(n));
         this->convertToTriplets(list, -2, ni, this->d_2(n));
         this->convertToTriplets(list, -1, ni, this->d_1(n));
         this->convertToTriplets(list, 0, ni, this->d0(n));
         this->convertToTriplets(list, 1, ni, this->d1(n));
         this->convertToTriplets(list, 2, ni, this->d2(n));
         this->convertToTriplets(list, 3, ni, this->d3(n));
         this->convertToTriplets(list, 4, ni, this->d4(n));
         this->convertToTriplets(list, 5, ni, this->d5(n));
         this->convertToTriplets(list, 6, ni, this->d6(n));
      }
   }

}
}
}
}
