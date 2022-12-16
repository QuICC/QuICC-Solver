/** 
 * @file R2.cpp
 * @brief Source of the implementation of the full sphere Worland R2 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Class include
//
#include "QuICC/SparseSM/Worland/R2.hpp"

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/R2Diags.hpp"
#include "QuICC/SparseSM/Worland/Legendre/R2Diags.hpp"
#include "QuICC/SparseSM/Worland/CylEnergy/R2Diags.hpp"
#include "QuICC/SparseSM/Worland/SphEnergy/R2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   R2::R2(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IWorlandOperator(rows, cols, alpha, dBeta)
   {
      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            this->mpImpl = std::make_shared<Chebyshev::R2Diags>(alpha, l, q);
            break;
         case WorlandKind::LEGENDRE:
            this->mpImpl = std::make_shared<Legendre::R2Diags>(alpha, l, q);
            break;
         case WorlandKind::CYLENERGY:
            this->mpImpl = std::make_shared<CylEnergy::R2Diags>(alpha, l, q);
            break;
         case WorlandKind::SPHENERGY:
            this->mpImpl = std::make_shared<SphEnergy::R2Diags>(alpha, l, q);
            break;
      }
   }

   R2::~R2()
   {
   }

   void R2::buildTriplets(ISparseSMOperator::TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows(), 0, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(5*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -1, ni.tail(ni.size()-1), this->mpImpl->d_1(n.tail(n.size()-1)));
         this->convertToTriplets(list, 0, ni, this->mpImpl->d0(n));
         this->convertToTriplets(list, 1, ni, this->mpImpl->d1(n));
      }
   }

}
}
}
