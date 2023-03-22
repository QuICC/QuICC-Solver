/** 
 * @file Id.cpp
 * @brief Source of the implementation of the Chebyshev Id sparse operator, y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Id.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   Id::Id(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const int q, const int s)
      : ILinearMapOperator(rows, cols, lower, upper), mQ(q), mShift(s)
   {
      int q_ = q;
      if(this->mShift < 0)
      {
         q_ = std::max(0, q + this->mShift);
      }

      this->mQ = q_;
   }

   Id::ACoeff_t Id::d0(const ACoeff_t& n) const
   {
      ACoeff_t val = ACoeff_t::Ones(n.size());

      if(this->mQ > 0)
      {
         val.topRows(this->mQ).setZero();
      }
      else if(this->mQ < 0)
      {
         val.bottomRows(-this->mQ).setZero();
      }

      return val;
   }

   void Id::buildTriplets(TripletList_t& list) const
   {
      int nN;
      if(this->mShift <= 0)
      {
         nN = this->rows();
      }
      else
      {
         nN = this->cols();
      }
      int n0 = std::max(0, -this->mShift);
      int nMax = nN-1;

      ACoeffI ni = ACoeffI::LinSpaced(nN - std::abs(this->mShift), n0, nMax);
      ACoeff_t n = (ni).cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, this->mShift, ni, this->d0(n));
      }
   }

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
