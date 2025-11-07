/**
 * @file Llm1DivS1.cpp
 * @brief Source of the implementation of the associated Legendre l(l+1)/Sin projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/Llm1DivS1.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   void Llm1DivS1<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      DivS1<base_t>::applyOperator(rOut, i, this->mLlm1.bottomRows(in.rows()).asDiagonal()*in);
   }

   void Llm1DivS1<base_t>::initSpecial() const
   {
      // Initialise storage for l(l+1) factor
      this->mLlm1 = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLlm1 = this->mLlm1.array()*(this->mLlm1.array() + 1.0)-1.0;
   }

}
}
}
}
}
