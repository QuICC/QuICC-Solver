/** 
 * @file I6CylLaplh.cpp
 * @brief Source of the implementation of the full sphere Worland I6CylLaplh sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/I6CylLaplh.hpp"

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I6CylLaplhDiags.hpp"
//#include "QuICC/SparseSM/Worland/Legendre/I6CylLaplhDiags.hpp"
//#include "QuICC/SparseSM/Worland/CylEnergy/I6CylLaplhDiags.hpp"
//#include "QuICC/SparseSM/Worland/SphEnergy/I6CylLaplhDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I6CylLaplh::I6CylLaplh(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IWorlandOperator(rows, cols, dBeta, l)
   {
      switch(this->type())
      {
         case CHEBYSHEV:
            this->mpImpl = std::make_shared<Chebyshev::I6CylLaplhDiags>(alpha, l);
            break;
         case LEGENDRE:
            throw std::logic_error("Operator is not implemented for Legendre type");
            //this->mpImpl = std::make_shared<Legendre::I6CylLaplhDiags>(alpha, l);
            break;
         case CYLENERGY:
            throw std::logic_error("Operator is not implemented for CylEnergy type");
            //this->mpImpl = std::make_shared<CylEnergy::I6CylLaplhDiags>(alpha, l);
            break;
         case SPHENERGY:
            throw std::logic_error("Operator is not implemented for SphEnergy type");
            //this->mpImpl = std::make_shared<SphEnergy::I6CylLaplhDiags>(alpha, l);
            break;
      }
   }

   I6CylLaplh::~I6CylLaplh()
   {
   }

   void I6CylLaplh::buildTriplets(TripletList_t& list) const
   {
      const int dShift = 3;
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-3, 3, this->rows()-1);
      ACoeff_t n = (ni + dShift).cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(11*std::max(this->rows(),this->cols()));
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
      }
   }

}
}
}
