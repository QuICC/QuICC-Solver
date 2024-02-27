/** 
 * @file InsulatingSphere.cpp
 * @brief Source of the implementation of the full sphere Worland insulating sphere sparse operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/InsulatingSphere.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Chebyshev/InsulatingSphereDiags.hpp"
//#include "QuICC/SparseSM/Worland/Stencil/Legendre/InsulatingSphereDiags.hpp"
//#include "QuICC/SparseSM/Worland/Stencil/CylEnergy/InsulatingSphereDiags.hpp"
#include "QuICC/SparseSM/Worland/Stencil/SphEnergy/InsulatingSphereDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

   InsulatingSphere::InsulatingSphere(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IWorlandOperator(rows, cols, alpha, dBeta)
   {
      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            this->mpImpl = std::make_shared<Chebyshev::InsulatingSphereDiags>(alpha, l);
            break;
         case WorlandKind::LEGENDRE:
            //this->mpImpl = std::make_shared<Legendre::InsulatingSphereDiags>(alpha, l);
            throw std::logic_error("Not yet implemented");
            break;
         case WorlandKind::CYLENERGY:
            //this->mpImpl = std::make_shared<CylEnergy::InsulatingSphereDiags>(alpha, l);
            throw std::logic_error("Not yet implemented");
            break;
         case WorlandKind::SPHENERGY:
            this->mpImpl = std::make_shared<SphEnergy::InsulatingSphereDiags>(alpha, l);
            break;
      }
   }

   void InsulatingSphere::buildTriplets(TripletList_t& list) const
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
