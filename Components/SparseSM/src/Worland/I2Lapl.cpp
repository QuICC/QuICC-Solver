/** 
 * @file I2Lapl.cpp
 * @brief Source of the implementation of the full sphere Worland I2Lapl sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/I2Lapl.hpp"

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I2LaplDiags.hpp"
//#include "QuICC/SparseSM/Worland/Legendre/I2LaplDiags.hpp"
//#include "QuICC/SparseSM/Worland/CylEnergy/I2LaplDiags.hpp"
//#include "QuICC/SparseSM/Worland/SphEnergy/I2LaplDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I2Lapl::I2Lapl(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IWorlandOperator(rows, cols, alpha, dBeta)
   {
      switch(this->type())
      {
         case CHEBYSHEV:
            this->mpImpl = std::make_shared<Chebyshev::I2LaplDiags>(alpha, l);
            break;
         case LEGENDRE:
            //this->mpImpl = std::make_shared<Legendre::I2LaplDiags>(alpha, l);
            throw std::logic_error("Not yet implemented");
            break;
         case CYLENERGY:
            //this->mpImpl = std::make_shared<CylEnergy::I2LaplDiags>(alpha, l);
            throw std::logic_error("Not yet implemented");
            break;
         case SPHENERGY:
            //this->mpImpl = std::make_shared<SphEnergy::I2LaplDiags>(alpha, l);
            throw std::logic_error("Not yet implemented");
            break;
      }
   }

   void I2Lapl::buildTriplets(TripletList_t& list) const
   {
      const int dShift = 1;
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-1, 1, this->rows()-1);
      ACoeff_t n = (ni + dShift).cast<Scalar_t>();

      if(n.size() > 0)
      {
         list.reserve(5*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -1 + dShift, ni, this->mpImpl->d_1(n));
         this->convertToTriplets(list, 0 + dShift, ni, this->mpImpl->d0(n));
         this->convertToTriplets(list, 1 + dShift, ni, this->mpImpl->d1(n));
      }
   }

   void I2Lapl::buildBanded(internal::Matrix& bd, unsigned int& kL, unsigned int &kU) const
   {
      kL = 1;
      kU = 3;
      bd.resize(kL+kU+1, this->rows());

      const int dShift = 1;
      ACoeffI ni = ACoeffI::LinSpaced(this->rows()-1, 1, this->rows()-1);
      ACoeff_t n = (ni + dShift).cast<Scalar_t>();

      int r = this->rows() - dShift;
      bd.row(1).rightCols(r-2) = this->mpImpl->d1(n).topRows(r-2);
      bd.block(1, 2, 1, dShift).setZero();
      bd.row(2).rightCols(r-1) = this->mpImpl->d0(n).topRows(r-1);
      bd.block(2, 1, 1, dShift).setZero();
      bd.row(3).rightCols(r) = this->mpImpl->d_1(n).topRows(r);
   }

}
}
}
