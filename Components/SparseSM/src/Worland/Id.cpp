/** 
 * @file Id.cpp
 * @brief Source of the implementation of the full sphere Worland Id sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Class include
//
#include "QuICC/SparseSM/Worland/Id.hpp"

// Project includes
//
#include "QuICC/SparseSM/Worland/IdDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   Id::Id(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q, const int s)
      : IWorlandOperator(rows, cols, alpha, dBeta), mShift(s)
   {
      int q_ = q;
      if(this->mShift < 0)
      {
         q_ = std::max(0, q + this->mShift);
      }
      this->mpImpl = std::make_shared<IdDiags>(alpha, dBeta, l, q_);
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
         this->convertToTriplets(list, this->mShift, ni, this->mpImpl->d0(n));
      }
   }

}
}
}
