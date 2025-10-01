/**
 * @file I3.cpp
 * @brief Source of the implementation of the full sphere Worland I3 sparse
 * operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I3Diags.hpp"
#include "QuICC/SparseSM/Worland/CylEnergy/I3Diags.hpp"
#include "QuICC/SparseSM/Worland/I3.hpp"
#include "QuICC/SparseSM/Worland/Legendre/I3Diags.hpp"
#include "QuICC/SparseSM/Worland/SphEnergy/I3Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

I3::I3(const int rows, const int cols, const Scalar_t alpha,
   const Scalar_t dBeta, const int l, const int q) :
    IWorlandOperator(rows, cols, alpha, dBeta)
{
   switch (this->type())
   {
   case WorlandKind::CHEBYSHEV:
      this->mpImpl = std::make_shared<Chebyshev::I3Diags>(alpha, l, q);
      break;
   case WorlandKind::LEGENDRE:
      this->mpImpl = std::make_shared<Legendre::I3Diags>(alpha, l, q);
      break;
   case WorlandKind::CYLENERGY:
      this->mpImpl = std::make_shared<CylEnergy::I3Diags>(alpha, l, q);
      break;
   case WorlandKind::SPHENERGY:
      this->mpImpl = std::make_shared<SphEnergy::I3Diags>(alpha, l, q);
      break;
   }
}

void I3::buildTriplets(TripletList_t& list) const
{
   const int dShift = 2;
   ACoeffI ni = ACoeffI::LinSpaced(this->rows() - 1, 1, this->rows() - 1);
   ACoeff_t n = (ni + dShift).cast<Scalar_t>();

   // Precompute the normalization factors
   int maxN = this->rows() - 1 + dShift + 3;
   this->mpImpl->precomputeNorm(maxN, 0);

   if (n.size() > 0)
   {
      list.reserve(7 * std::max(this->rows(), this->cols()));
      this->convertToTriplets(list, -3 + dShift, ni, this->mpImpl->d_3(n));
      this->convertToTriplets(list, -2 + dShift, ni, this->mpImpl->d_2(n));
      this->convertToTriplets(list, -1 + dShift, ni, this->mpImpl->d_1(n));
      this->convertToTriplets(list, 0 + dShift, ni, this->mpImpl->d0(n));
      this->convertToTriplets(list, 1 + dShift, ni, this->mpImpl->d1(n));
      this->convertToTriplets(list, 2 + dShift, ni, this->mpImpl->d2(n));
      this->convertToTriplets(list, 3 + dShift, ni, this->mpImpl->d3(n));
   }
}

void I3::buildBanded(Internal::Matrix& bd, unsigned int& kL,
   unsigned int& kU) const
{
   kL = 1;
   kU = 5;
   bd.resize(kL + kU + 1, this->rows());

   const int dShift = 2;
   ACoeffI ni = ACoeffI::LinSpaced(this->rows() - 2, 2, this->rows() - 1);
   ACoeff_t n = (ni + dShift).cast<Scalar_t>();

   int r = this->rows() - dShift;
   bd.row(1).rightCols(r - 5) = this->mpImpl->d3(n).topRows(r - 5);
   bd.block(1, 5, 1, dShift).setZero();
   bd.row(2).rightCols(r - 4) = this->mpImpl->d2(n).topRows(r - 4);
   bd.block(2, 4, 1, dShift).setZero();
   bd.row(3).rightCols(r - 3) = this->mpImpl->d1(n).topRows(r - 3);
   bd.block(3, 3, 1, dShift).setZero();
   bd.row(4).rightCols(r - 2) = this->mpImpl->d0(n).topRows(r - 2);
   bd.block(4, 2, 1, dShift).setZero();
   bd.row(5).rightCols(r - 1) = this->mpImpl->d_1(n).topRows(r - 1);
   bd.block(5, 1, 1, dShift).setZero();
   bd.row(6).rightCols(r) = this->mpImpl->d_2(n).topRows(r);
   bd.block(6, 0, 1, dShift).setZero();
   bd.block(7, 1, 1, r) = this->mpImpl->d_3(n).topRows(r).transpose();
}

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
