/** 
 * @file Id.cpp
 * @brief Source of the implementation of the full sphere Worland Id sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/Id.hpp"

// Project includes
//
#include "QuICC/SparseSM/Worland/IdDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   Id::Id(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IWorlandOperator(rows, cols, alpha, dBeta)
   {
      this->mpImpl = std::make_shared<IdDiags>(alpha, dBeta, l, q);
   }

   void Id::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows(), 0, this->rows()-1);
      ACoeff_t n = (ni).cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, 0, ni, this->mpImpl->d0(n));
      }
   }

}
}
}
