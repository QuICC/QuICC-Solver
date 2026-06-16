/**
 * @file I3D3p.cpp
 * @brief Source of the implementation of the full sphere Worland I3D3p sparse
 * operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I3D3pDiags.hpp"
#include "QuICC/SparseSM/Worland/I3D3p.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

I3D3p::I3D3p(const int rows, const int cols, const Scalar_t alpha,
   const Scalar_t dBeta, const int l, const int q) :
    IWorlandOperator(rows, cols, alpha, dBeta)
{
   switch (this->type())
   {
   case WorlandKind::CHEBYSHEV:
      this->mpImpl = std::make_shared<Chebyshev::I3D3pDiags>(alpha, l, q);
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

void I3D3p::buildTriplets(TripletList_t& list) const
{
   const int dShift = 2;
   ACoeffI ni = ACoeffI::LinSpaced(this->rows() - 1, 1, this->rows() - 1);
   ACoeff_t n = (ni + dShift).cast<Scalar_t>();

   // Precompute the normalization factors (base l and the l+3 coupling)
   int maxN = this->rows() - 1 + dShift + 3;
   this->mpImpl->precomputeNorm(maxN, 0);
   this->mpImpl->precomputeNorm(maxN, 3);

   if (n.size() > 0)
   {
      list.reserve(4 * std::max(this->rows(), this->cols()));
      this->convertToTriplets(list, -3 + dShift, ni, this->mpImpl->d_3(n));
      this->convertToTriplets(list, -2 + dShift, ni, this->mpImpl->d_2(n));
      this->convertToTriplets(list, -1 + dShift, ni, this->mpImpl->d_1(n));
      this->convertToTriplets(list, 0 + dShift, ni, this->mpImpl->d0(n));
   }
}

void I3D3p::buildBanded(Internal::Matrix& bd, unsigned int& kL,
   unsigned int& kU) const
{
   throw std::logic_error("Banded matrix is not yet implemented");
}

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
