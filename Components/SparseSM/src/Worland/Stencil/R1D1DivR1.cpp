/** 
 * @file R1D1DivR1.cpp
 * @brief Source of the implementation of the full sphere Worland R1D1DivR1 sparse operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Chebyshev/R1D1DivR1Diags.hpp"
//#include "QuICC/SparseSM/Worland/Stencil/Legendre/R1D1DivR1Diags.hpp"
//#include "QuICC/SparseSM/Worland/Stencil/CylEnergy/R1D1DivR1Diags.hpp"
#include "QuICC/SparseSM/Worland/Stencil/SphEnergy/R1D1DivR1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

   R1D1DivR1::R1D1DivR1(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IWorlandOperator(rows, cols, alpha, dBeta)
   {
      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            this->mpImpl = std::make_shared<Chebyshev::R1D1DivR1Diags>(alpha, l);
            break;
         case WorlandKind::LEGENDRE:
            //this->mpImpl = std::make_shared<Legendre::R1D1DivR1Diags>(alpha, l);
            throw std::logic_error("Not yet implemented");
            break;
         case WorlandKind::CYLENERGY:
            //this->mpImpl = std::make_shared<CylEnergy::R1D1DivR1Diags>(alpha, l);
            throw std::logic_error("Not yet implemented");
            break;
         case WorlandKind::SPHENERGY:
            this->mpImpl = std::make_shared<SphEnergy::R1D1DivR1Diags>(alpha, l);
            break;
      }
   }

   void R1D1DivR1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows(), 0, this->rows()-1);
      ACoeff_t n = (ni).cast<Scalar_t>();
      ACoeffI ni1 = ni.bottomRows(this->rows()-1);
      ACoeff_t n1 = n.bottomRows(this->rows()-1);

      if(n.size() > 0)
      {
         list.reserve(2*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -1, ni1, this->mpImpl->d_1(n1));
         this->convertToTriplets(list, 0, ni, this->mpImpl->d0(n));
      }
   }

} // Stencil
} // Worland
} // SparsesM
} // QuICC
