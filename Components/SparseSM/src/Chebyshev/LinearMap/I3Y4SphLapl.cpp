/**
 * @file I3Y4SphLapl.cpp
 * @brief Source of the implementation of the I^3 Y^4 spherical laplacian sparse operator, with y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I3Y4SphLapl.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   using namespace Internal::Literals;

   I3Y4SphLapl::I3Y4SphLapl(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t l)
      : ISphericalOperator(rows, cols, lower, upper, l)
   {
   }

   I3Y4SphLapl::ACoeff_t I3Y4SphLapl::d_5(const ACoeff_t& n) const
   {
      return -a<5>()*(l<1>() - n + 5.0_mp)*(l<1>() + n - 4.0_mp)/(32.0_mp*n*(n - 2.0_mp)*(n - 1.0_mp));
   }

   I3Y4SphLapl::ACoeff_t I3Y4SphLapl::d_4(const ACoeff_t& n) const
   {
      return -a<4>()*b<1>()*(l<2>() + l<1>() - 2.0_mp*n*n + 15.0_mp*n - 28.0_mp)/(8.0_mp*n*(n - 2.0_mp)*(n - 1.0_mp));
   }

   I3Y4SphLapl::ACoeff_t I3Y4SphLapl::d_3(const ACoeff_t& n) const
   {
      return a<3>()*(a<2>()*l<2>()*n - 5.0_mp*a<2>()*l<2>() + a<2>()*l<1>()*n - 5.0_mp*a<2>()*l<1>() + 3.0_mp*a<2>()*n*n*n - 12.0_mp*a<2>()*n*n - 15.0_mp*a<2>()*n + 72.0_mp*a<2>() - 4.0_mp*b<2>()*l<2>()*n - 4.0_mp*b<2>()*l<2>() - 4.0_mp*b<2>()*l<1>()*n - 4.0_mp*b<2>()*l<1>() + 24.0_mp*b<2>()*n*n*n - 120.0_mp*b<2>()*n*n + 72.0_mp*b<2>()*n + 216.0_mp*b<2>())/(32.0_mp*n*(n - 2.0_mp)*(n - 1.0_mp)*(n + 1.0_mp));
   }

   I3Y4SphLapl::ACoeff_t I3Y4SphLapl::d_2(const ACoeff_t& n) const
   {
      return a<2>()*b<1>()*(a<2>()*l<2>() + a<2>()*l<1>() + 2.0_mp*a<2>()*n*n - 17.0_mp*a<2>() + 4.0_mp*b<2>()*n*n - 6.0_mp*b<2>()*n - 10.0_mp*b<2>())/(4.0_mp*n*(n - 1.0_mp)*(n + 1.0_mp));
   }

   I3Y4SphLapl::ACoeff_t I3Y4SphLapl::d_1(const ACoeff_t& n) const
   {
      return a<1>()*(a<4>()*l<2>()*n + 4.0_mp*a<4>()*l<2>() + a<4>()*l<1>()*n + 4.0_mp*a<4>()*l<1>() + a<4>()*n*n*n + 7.0_mp*a<4>()*n*n - 10.0_mp*a<4>()*n - 52.0_mp*a<4>() + 6.0_mp*a<2>()*b<2>()*l<2>()*n + 12.0_mp*a<2>()*b<2>()*l<2>() + 6.0_mp*a<2>()*b<2>()*l<1>()*n + 12.0_mp*a<2>()*b<2>()*l<1>() + 12.0_mp*a<2>()*b<2>()*n*n*n + 48.0_mp*a<2>()*b<2>()*n*n - 84.0_mp*a<2>()*b<2>()*n - 264.0_mp*a<2>()*b<2>() + 8.0_mp*b<4>()*n*n*n + 8.0_mp*b<4>()*n*n - 32.0_mp*b<4>()*n - 32.0_mp*b<4>())/(16.0_mp*n*(n - 2.0_mp)*(n + 1.0_mp)*(n + 2.0_mp));
   }

   I3Y4SphLapl::ACoeff_t I3Y4SphLapl::d0(const ACoeff_t& n) const
   {
      return 3.0_mp*a<2>()*b<1>()*(a<2>()*l<2>() + a<2>()*l<1>() + 3.0_mp*a<2>()*n*n - 18.0_mp*a<2>() + 4.0_mp*b<2>()*n*n - 16.0_mp*b<2>())/(4.0_mp*(n - 2.0_mp)*(n - 1.0_mp)*(n + 1.0_mp)*(n + 2.0_mp));
   }

   I3Y4SphLapl::ACoeff_t I3Y4SphLapl::d1(const ACoeff_t& n) const
   {
      return -a<1>()*(a<4>()*l<2>()*n - 4.0_mp*a<4>()*l<2>() + a<4>()*l<1>()*n - 4.0_mp*a<4>()*l<1>() + a<4>()*n*n*n - 7.0_mp*a<4>()*n*n - 10.0_mp*a<4>()*n + 52.0_mp*a<4>() + 6.0_mp*a<2>()*b<2>()*l<2>()*n - 12.0_mp*a<2>()*b<2>()*l<2>() + 6.0_mp*a<2>()*b<2>()*l<1>()*n - 12.0_mp*a<2>()*b<2>()*l<1>() + 12.0_mp*a<2>()*b<2>()*n*n*n - 48.0_mp*a<2>()*b<2>()*n*n - 84.0_mp*a<2>()*b<2>()*n + 264.0_mp*a<2>()*b<2>() + 8.0_mp*b<4>()*n*n*n - 8.0_mp*b<4>()*n*n - 32.0_mp*b<4>()*n + 32.0_mp*b<4>())/(16.0_mp*n*(n - 2.0_mp)*(n - 1.0_mp)*(n + 2.0_mp));
   }

   I3Y4SphLapl::ACoeff_t I3Y4SphLapl::d2(const ACoeff_t& n) const
   {
      return -a<2>()*b<1>()*(a<2>()*l<2>() + a<2>()*l<1>() + 2.0_mp*a<2>()*n*n - 17.0_mp*a<2>() + 4.0_mp*b<2>()*n*n + 6.0_mp*b<2>()*n - 10.0_mp*b<2>())/(4.0_mp*n*(n - 1.0_mp)*(n + 1.0_mp));
   }

   I3Y4SphLapl::ACoeff_t I3Y4SphLapl::d3(const ACoeff_t& n) const
   {
      return -a<3>()*(a<2>()*l<2>()*n + 5.0_mp*a<2>()*l<2>() + a<2>()*l<1>()*n + 5.0_mp*a<2>()*l<1>() + 3.0_mp*a<2>()*n*n*n + 12.0_mp*a<2>()*n*n - 15.0_mp*a<2>()*n - 72.0_mp*a<2>() - 4.0_mp*b<2>()*l<2>()*n + 4.0_mp*b<2>()*l<2>() - 4.0_mp*b<2>()*l<1>()*n + 4.0_mp*b<2>()*l<1>() + 24.0_mp*b<2>()*n*n*n + 120.0_mp*b<2>()*n*n + 72.0_mp*b<2>()*n - 216.0_mp*b<2>())/(32.0_mp*n*(n - 1.0_mp)*(n + 1.0_mp)*(n + 2.0_mp));
   }

   I3Y4SphLapl::ACoeff_t I3Y4SphLapl::d4(const ACoeff_t& n) const
   {
      return a<4>()*b<1>()*(l<2>() + l<1>() - 2.0_mp*n*n - 15.0_mp*n - 28.0_mp)/(8.0_mp*n*(n + 1.0_mp)*(n + 2.0_mp));
   }

   I3Y4SphLapl::ACoeff_t I3Y4SphLapl::d5(const ACoeff_t& n) const
   {
      return a<5>()*(l<1>() - n - 4.0_mp)*(l<1>() + n + 5.0_mp)/(32.0_mp*n*(n + 1.0_mp)*(n + 2.0_mp));
   }

   void I3Y4SphLapl::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-3, 3, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(11*std::max(this->rows(),this->cols()));
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
      }
   }

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
