/**
 * @file Id.cpp
 * @brief Source of the implementation of the full sphere Bessel Id sparse
 * operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/Id.hpp"
#include "QuICC/SparseSM/Bessel/IdDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

Id::Id(const int rows, const int cols, const BesselKind type, const int l,
   const int s) :
    IBesselOperator(rows, cols, type), mShift(s)
{
   this->mpImpl = std::make_shared<IdDiags>(type, l);
}

void Id::buildTriplets(TripletList_t& list) const
{
   int nN;
   if (this->mShift <= 0)
   {
      nN = this->rows();
   }
   else
   {
      nN = this->cols();
   }
   int n0 = std::max(0, -this->mShift);
   int nMax = nN - 1;

   ACoeffI ni = ACoeffI::LinSpaced(nN - std::abs(this->mShift), n0, nMax);
   ACoeff_t n = (ni).cast<Scalar_t>();

   if (n.size() > 0)
   {
      list.reserve(std::max(this->rows(), this->cols()));
      this->convertToTriplets(list, this->mShift, ni, this->mpImpl->d0(n));
   }
}

} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC
