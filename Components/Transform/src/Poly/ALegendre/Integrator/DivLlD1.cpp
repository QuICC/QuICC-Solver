/**
 * @file DivLlD1.cpp
 * @brief Source of the implementation of the associated Legendre 1/l(l+1) D integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/DivLlD1.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   DivLlD1::DivLlD1()
   {
   }

   DivLlD1::~DivLlD1()
   {
   }

   void DivLlD1::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      D1::applyOperator(rOut, i, in);
      rOut = this->mDivLl.bottomRows(rOut.rows()).asDiagonal()*rOut;
   }

   void DivLlD1::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mDivLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mDivLl = (this->mDivLl.array()*(this->mDivLl.array() + 1.0)).pow(-1);
      this->mDivLl(0) = 0.0;
   }

}
}
}
}
}
