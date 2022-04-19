/**
 * @file DivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre 1/S D_phi projector
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Projector/DivS1Dp.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   DivS1Dp::DivS1Dp()
   {
   }

   DivS1Dp::~DivS1Dp()
   {
   }

   void DivS1Dp::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      MHDComplex c(0.0,static_cast<MHDFloat>(this->mspSetup->slow(i)));
      DivS1::applyOperator(rOut, i, c*in);
   }

}
}
}
}
}
