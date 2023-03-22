/** 
 * @file ValueD1.cpp
 * @brief Source of the implementation of the full sphere Worland ValueD1 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/ValueD1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Chebyshev/ValueD1Diags.hpp"
//#include "QuICC/SparseSM/Worland/Legendre/ValueD1Diags.hpp"
//#include "QuICC/SparseSM/Worland/CylEnergy/ValueD1Diags.hpp"
//#include "QuICC/SparseSM/Worland/SphEnergy/ValueD1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

   ValueD1::ValueD1(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IWorlandOperator(rows, cols, alpha, dBeta)
   {
      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            this->mpImpl = std::make_shared<Chebyshev::ValueD1Diags>(alpha, l);
            break;
         case WorlandKind::LEGENDRE:
            //this->mpImpl = std::make_shared<Legendre::ValueD1Diags>(alpha, l);
            throw std::logic_error("Not yet implemented");
            break;
         case WorlandKind::CYLENERGY:
            //this->mpImpl = std::make_shared<CylEnergy::ValueD1Diags>(alpha, l);
            throw std::logic_error("Not yet implemented");
            break;
         case WorlandKind::SPHENERGY:
            //this->mpImpl = std::make_shared<SphEnergy::ValueD1Diags>(alpha, l);
            throw std::logic_error("Not yet implemented");
            break;
      }
   }

   void ValueD1::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows(), 0, this->rows()-1);
      ACoeff_t n = (ni).cast<Scalar_t>();
      ACoeffI ni1 = ni.bottomRows(this->rows()-1);
      ACoeff_t n1 = n.bottomRows(this->rows()-1);
      ACoeffI ni2 = ni.bottomRows(this->rows()-2);
      ACoeff_t n2 = n.bottomRows(this->rows()-2);

      if(n.size() > 0)
      {
         list.reserve(3*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -2, ni2, this->mpImpl->d_2(n2));
         this->convertToTriplets(list, -1, ni1, this->mpImpl->d_1(n1));
         this->convertToTriplets(list, 0, ni, this->mpImpl->d0(n));
      }
   }

} // Stencil
} // Worland
} // SparsesM
} // QuICC
