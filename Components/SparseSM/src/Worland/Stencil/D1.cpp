/** 
 * @file D1.cpp
 * @brief Source of the implementation of the full sphere Worland D1 sparse operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/D1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Chebyshev/D1Diags.hpp"
//#include "QuICC/SparseSM/Worland/Stencil/Legendre/D1Diags.hpp"
//#include "QuICC/SparseSM/Worland/Stencil/CylEnergy/D1Diags.hpp"
#include "QuICC/SparseSM/Worland/Stencil/SphEnergy/D1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

   D1::D1(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IWorlandOperator(rows, cols, alpha, dBeta)
   {
      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            this->mpImpl = std::make_shared<Chebyshev::D1Diags>(alpha, l);
            break;
         case WorlandKind::LEGENDRE:
            //this->mpImpl = std::make_shared<Legendre::D1Diags>(alpha, l);
            throw std::logic_error("Not yet implemented");
            break;
         case WorlandKind::CYLENERGY:
            //this->mpImpl = std::make_shared<CylEnergy::D1Diags>(alpha, l);
            throw std::logic_error("Not yet implemented");
            break;
         case WorlandKind::SPHENERGY:
            this->mpImpl = std::make_shared<SphEnergy::D1Diags>(alpha, l);
            break;
      }
   }

   void D1::buildTriplets(TripletList_t& list) const
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
