/**
 * @file LlDivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1) 1/Sin D_phi integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/LlDivS1Dp.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   LlDivS1Dp::LlDivS1Dp()
   {
   }

   LlDivS1Dp::~LlDivS1Dp()
   {
   }

   void LlDivS1Dp::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      DivS1::applyOperator(rOut, i, in);
      rOut = (static_cast<MHDFloat>(-this->mspSetup->slow(i))*Math::cI*this->mLl.bottomRows(rOut.rows())).asDiagonal()*rOut;
   }

   void LlDivS1Dp::initSpecial() const
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
