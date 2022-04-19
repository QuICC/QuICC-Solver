/**
 * @file Ll2.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1)^2 P integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Ll2.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   Ll2::Ll2()
   {
   }

   Ll2::~Ll2()
   {
   }

   void Ll2::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      P::applyOperator(rOut, i, in);
      rOut = this->mLl2.bottomRows(rOut.rows()).asDiagonal()*rOut;
   }

   void Ll2::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLl2 = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl2 = (this->mLl2.array()*(this->mLl2.array() + 1.0)).pow(2);
   }

}
}
}
}
}
