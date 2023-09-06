/** 
 * @file I4D1R1.cpp
 * @brief Source of the implementation of the full sphere Worland I4D1R1 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I4D1R1.hpp"
#include "QuICC/SparseSM/Worland/Chebyshev/I4D1R1Diags.hpp"
//#include "QuICC/SparseSM/Worland/Legendre/I4D1R1Diags.hpp"
//#include "QuICC/SparseSM/Worland/CylEnergy/I4D1R1Diags.hpp"
//#include "QuICC/SparseSM/Worland/SphEnergy/I4D1R1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I4D1R1::I4D1R1(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IWorlandOperator(rows, cols, alpha, dBeta)
   {
      if(l != 1)
      {
         throw std::logic_error("Operator only exists for l = 1");
      }

      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            this->mpImpl = std::make_shared<Chebyshev::I4D1R1Diags>(alpha, l, q);
            break;
         case WorlandKind::LEGENDRE:
            throw std::logic_error("Operator is not implemented for Legendre type");
            //this->mpImpl = std::make_shared<Legendre::I4D1R1Diags>(alpha, l, q);
            break;
         case WorlandKind::CYLENERGY:
            throw std::logic_error("Operator is not implemented for CylEnergy type");
            //this->mpImpl = std::make_shared<CylEnergy::I4D1R1Diags>(alpha, l, q);
            break;
         case WorlandKind::SPHENERGY:
            throw std::logic_error("Operator is not implemented for SphEnergy type");
            //this->mpImpl = std::make_shared<SphEnergy::I4D1R1Diags>(alpha, l, q);
            break;
      }
   }

   void I4D1R1::buildTriplets(TripletList_t& list) const
   {
      const int dShift = 2;
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-2, 2, this->rows()-1);
      ACoeff_t n = (ni + dShift).cast<Scalar_t>();

      // Precompute the normalization factors
      int maxN = this->rows()-1 + dShift + 4;
      this->mpImpl->precomputeNorm(maxN, 0);

      if(n.size() > 0)
      {
         list.reserve(8*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -3 + dShift, ni, this->mpImpl->d_3(n));
         this->convertToTriplets(list, -2 + dShift, ni, this->mpImpl->d_2(n));
         this->convertToTriplets(list, -1 + dShift, ni, this->mpImpl->d_1(n));
         this->convertToTriplets(list, 0 + dShift, ni, this->mpImpl->d0(n));
         this->convertToTriplets(list, 1 + dShift, ni, this->mpImpl->d1(n));
         this->convertToTriplets(list, 2 + dShift, ni, this->mpImpl->d2(n));
         this->convertToTriplets(list, 3 + dShift, ni, this->mpImpl->d3(n));
         this->convertToTriplets(list, 4 + dShift, ni, this->mpImpl->d4(n));
      }
   }

}
}
}
