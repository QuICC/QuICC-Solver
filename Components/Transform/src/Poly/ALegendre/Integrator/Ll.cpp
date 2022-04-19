/**
 * @file Ll.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1) P integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Ll.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   Ll::Ll()
   {
   }

   Ll::~Ll()
   {
   }

   void Ll::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      P::applyOperator(rOut, i, in);
      rOut = this->mLl.bottomRows(rOut.rows()).asDiagonal()*rOut;
   }

   void Ll::initSpecial() const
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
