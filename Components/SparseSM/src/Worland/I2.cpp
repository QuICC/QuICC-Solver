/** 
 * @file I2.cpp
 * @brief Source of the implementation of the full sphere Worland I2 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/WorlandKind.hpp"
#include "QuICC/SparseSM/Worland/Chebyshev/I2Diags.hpp"
#include "QuICC/SparseSM/Worland/Legendre/I2Diags.hpp"
#include "QuICC/SparseSM/Worland/CylEnergy/I2Diags.hpp"
#include "QuICC/SparseSM/Worland/SphEnergy/I2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I2::I2(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IWorlandOperator(rows, cols, alpha, dBeta)
   {
      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            this->mpImpl = std::make_shared<Chebyshev::I2Diags>(alpha, l, q);
            break;
         case WorlandKind::LEGENDRE:
            this->mpImpl = std::make_shared<Legendre::I2Diags>(alpha, l, q);
            break;
         case WorlandKind::CYLENERGY:
            this->mpImpl = std::make_shared<CylEnergy::I2Diags>(alpha, l, q);
            break;
         case WorlandKind::SPHENERGY:
            this->mpImpl = std::make_shared<SphEnergy::I2Diags>(alpha, l, q);
            break;
      }
   }

   void I2::buildTriplets(TripletList_t& list) const
   {
      const int dShift = 1;
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-1, 1, this->rows()-1);
      ACoeff_t n = (ni + dShift).cast<Scalar_t>();

      // Precompute the normalization factors
      int maxN = this->rows()-1 + dShift + 2;
      this->mpImpl->precomputeNorm(maxN, 0);

      if(n.size() > 0)
      {
         list.reserve(5*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -2 + dShift, ni, this->mpImpl->d_2(n));
         this->convertToTriplets(list, -1 + dShift, ni, this->mpImpl->d_1(n));
         this->convertToTriplets(list, 0 + dShift, ni, this->mpImpl->d0(n));
         this->convertToTriplets(list, 1 + dShift, ni, this->mpImpl->d1(n));
         this->convertToTriplets(list, 2 + dShift, ni, this->mpImpl->d2(n));
      }
   }

   void I2::buildBanded(internal::Matrix& bd, unsigned int& kL, unsigned int &kU) const
   {
      kL = 1;
      kU = 3;
      bd.resize(kL+kU+1, this->rows());

      const int dShift = 1;
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-1, 1, this->rows()-1);
      ACoeff_t n = (ni + dShift).cast<Scalar_t>();

      int r = this->rows() - dShift;
      bd.row(0).rightCols(r-3) = this->mpImpl->d2(n).topRows(r-3);
      bd.block(0, 3, 1, dShift).setZero();
      bd.row(1).rightCols(r-2) = this->mpImpl->d1(n).topRows(r-2);
      bd.block(1, 2, 1, dShift).setZero();
      bd.row(2).rightCols(r-1) = this->mpImpl->d0(n).topRows(r-1);
      bd.block(2, 1, 1, dShift).setZero();
      bd.row(3).rightCols(r) = this->mpImpl->d_1(n).topRows(r);
      bd.block(3, 0, 1, dShift).setZero();
      bd.row(4).leftCols(r) = this->mpImpl->d_2(n).topRows(r);
   }

}
}
}
