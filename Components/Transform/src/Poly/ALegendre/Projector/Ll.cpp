/**
 * @file Ll.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1) P projector
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Ll.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   Ll::Ll()
   {
   }

   Ll::~Ll()
   {
   }

   void Ll::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      P::applyOperator(rOut, i, this->mLl.bottomRows(in.rows()).asDiagonal()*in);
   }

   void Ll::initSpecial() const
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
