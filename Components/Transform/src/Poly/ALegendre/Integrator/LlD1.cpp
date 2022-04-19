/**
 * @file LlD1.cpp
 * @brief Source of the implementation of the associated Legendre 1/l(l+1) D integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/LlD1.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   LlD1::LlD1()
   {
   }

   LlD1::~LlD1()
   {
   }

   void LlD1::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      D1::applyOperator(rOut, i, in);
      rOut = this->mLl.bottomRows(rOut.rows()).asDiagonal()*rOut;
   }

   void LlD1::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl = (this->mLl.array()*(this->mLl.array() + 1.0));
   }

}
}
}
}
}
