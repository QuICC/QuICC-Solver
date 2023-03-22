/** 
 * @file I6.cpp
 * @brief Source of the implementation of the full sphere Worland I6 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I6.hpp"
#include "QuICC/SparseSM/Worland/Chebyshev/I6Diags.hpp"
#include "QuICC/SparseSM/Worland/Legendre/I6Diags.hpp"
#include "QuICC/SparseSM/Worland/CylEnergy/I6Diags.hpp"
#include "QuICC/SparseSM/Worland/SphEnergy/I6Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I6::I6(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IWorlandOperator(rows, cols, alpha, dBeta)
   {
      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            this->mpImpl = std::make_shared<Chebyshev::I6Diags>(alpha, l, q);
            break;
         case WorlandKind::LEGENDRE:
            this->mpImpl = std::make_shared<Legendre::I6Diags>(alpha, l, q);
            break;
         case WorlandKind::CYLENERGY:
            this->mpImpl = std::make_shared<CylEnergy::I6Diags>(alpha, l, q);
            break;
         case WorlandKind::SPHENERGY:
            this->mpImpl = std::make_shared<SphEnergy::I6Diags>(alpha, l, q);
            break;
      }
   }

   void I6::buildTriplets(TripletList_t& list) const
   {
      const int dShift = 3;
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-3, 3, this->rows()-1);
      ACoeff_t n = (ni + dShift).cast<Scalar_t>();

      // Precompute the normalization factors
      int maxN = this->rows()-1 + dShift + 6;
      this->mpImpl->precomputeNorm(maxN, 0);

      if(n.size() > 0)
      {
         list.reserve(13*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -6 + dShift, ni, this->mpImpl->d_6(n));
         this->convertToTriplets(list, -5 + dShift, ni, this->mpImpl->d_5(n));
         this->convertToTriplets(list, -4 + dShift, ni, this->mpImpl->d_4(n));
         this->convertToTriplets(list, -3 + dShift, ni, this->mpImpl->d_3(n));
         this->convertToTriplets(list, -2 + dShift, ni, this->mpImpl->d_2(n));
         this->convertToTriplets(list, -1 + dShift, ni, this->mpImpl->d_1(n));
         this->convertToTriplets(list, 0 + dShift, ni, this->mpImpl->d0(n));
         this->convertToTriplets(list, 1 + dShift, ni, this->mpImpl->d1(n));
         this->convertToTriplets(list, 2 + dShift, ni, this->mpImpl->d2(n));
         this->convertToTriplets(list, 3 + dShift, ni, this->mpImpl->d3(n));
         this->convertToTriplets(list, 4 + dShift, ni, this->mpImpl->d4(n));
         this->convertToTriplets(list, 5 + dShift, ni, this->mpImpl->d5(n));
         this->convertToTriplets(list, 6 + dShift, ni, this->mpImpl->d6(n));
      }
   }

}
}
}
