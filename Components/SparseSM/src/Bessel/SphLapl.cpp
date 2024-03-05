/**
 * @file SphLapl.cpp
 * @brief Source of the implementation of the full sphere Bessel SphLapl sparse
 * operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/SphLapl.hpp"
#include "QuICC/SparseSM/Bessel/Insulating/SphLaplDiags.hpp"
#include "QuICC/SparseSM/Bessel/Value/SphLaplDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

SphLapl::SphLapl(const int rows, const int cols, const BesselKind type,
   const int l) :
    IBesselOperator(rows, cols, type)
{
   switch (this->type())
   {
   case BesselKind::VALUE:
      this->mpImpl = std::make_shared<Value::SphLaplDiags>(l);
      break;
   case BesselKind::INSULATING:
      this->mpImpl = std::make_shared<Insulating::SphLaplDiags>(l);
      break;
   }
}

void SphLapl::buildTriplets(TripletList_t& list) const
{
   ACoeffI ni = ACoeffI::LinSpaced(this->rows(), 0, this->rows() - 1);
   ACoeff_t n = ni.cast<Scalar_t>();

   if (n.size() > 0)
   {
      list.reserve(1 * std::max(this->rows(), this->cols()));
      this->convertToTriplets(list, 0, ni, this->mpImpl->d0(n));
   }
}

void SphLapl::buildBanded(Internal::Matrix& bd, unsigned int& kL,
   unsigned int& kU) const
{
   kL = 1;
   kU = 3;
   bd.resize(kL + kU + 1, this->rows());

   const int dShift = 0;
   ACoeffI ni = ACoeffI::LinSpaced(this->rows() - 1, 1, this->rows() - 1);
   ACoeff_t n = (ni + dShift).cast<Scalar_t>();

   int r = this->rows() - dShift;
   bd.row(2).rightCols(r - 1) = this->mpImpl->d0(n).topRows(r - 1);
   bd.block(2, 1, 1, dShift).setZero();
}

} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC
