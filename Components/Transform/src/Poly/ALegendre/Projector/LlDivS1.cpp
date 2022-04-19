/**
 * @file LlDivS1.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1)/Sin projector
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Projector/LlDivS1.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   LlDivS1::LlDivS1()
   {
   }

   LlDivS1::~LlDivS1()
   {
   }

   void LlDivS1::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      DivS1::applyOperator(rOut, i, this->mLl.bottomRows(in.rows()).asDiagonal()*in);
   }

   void LlDivS1::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl = this->mLl.array()*(this->mLl.array() + 1.0);
   }

}
}
}
}
}
