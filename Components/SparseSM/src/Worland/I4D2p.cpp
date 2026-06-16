/**
 * @file I4D2p.cpp
 * @brief Source of the implementation of the full sphere Worland I4D2p sparse
 * operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I4D2pDiags.hpp"
#include "QuICC/SparseSM/Worland/I4D2p.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

I4D2p::I4D2p(const int rows, const int cols, const Scalar_t alpha,
   const Scalar_t dBeta, const int l, const int q) :
    IWorlandOperator(rows, cols, alpha, dBeta)
{
   switch (this->type())
   {
   case WorlandKind::CHEBYSHEV:
      this->mpImpl = std::make_shared<Chebyshev::I4D2pDiags>(alpha, l, q);
      break;
   case WorlandKind::LEGENDRE:
      throw std::logic_error("Not yet implemented");
      break;
   case WorlandKind::CYLENERGY:
      throw std::logic_error("Not yet implemented");
      break;
   case WorlandKind::SPHENERGY:
      throw std::logic_error("Not yet implemented");
      break;
   }
}

void I4D2p::buildTriplets(TripletList_t& list) const
{
   const int dShift = 2;
   ACoeffI ni = ACoeffI::LinSpaced(this->rows() - 2, 2, this->rows() - 1);
   ACoeff_t n = (ni + dShift).cast<Scalar_t>();

   // Precompute the normalization factors (base l and the l+2 coupling)
   int maxN = this->rows() - 1 + dShift + 5;
   this->mpImpl->precomputeNorm(maxN, 0);
   this->mpImpl->precomputeNorm(maxN, 2);

   if (n.size() > 0)
   {
      list.reserve(7 * std::max(this->rows(), this->cols()));
      this->convertToTriplets(list, -4 + dShift, ni, this->mpImpl->d_4(n));
      this->convertToTriplets(list, -3 + dShift, ni, this->mpImpl->d_3(n));
      this->convertToTriplets(list, -2 + dShift, ni, this->mpImpl->d_2(n));
      this->convertToTriplets(list, -1 + dShift, ni, this->mpImpl->d_1(n));
      this->convertToTriplets(list, 0 + dShift, ni, this->mpImpl->d0(n));
      this->convertToTriplets(list, 1 + dShift, ni, this->mpImpl->d1(n));
      this->convertToTriplets(list, 2 + dShift, ni, this->mpImpl->d2(n));
   }
}

void I4D2p::buildBanded(Internal::Matrix& bd, unsigned int& kL,
   unsigned int& kU) const
{
   throw std::logic_error("Banded matrix is not yet implemented");
}

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
