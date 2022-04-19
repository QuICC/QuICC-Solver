/**
 * @file DivS1D1S1.cpp
 * @brief Source of the implementation of the associated Legendre 1/Sin D Sin projector
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Projector/DivS1D1S1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   DivS1D1S1::DivS1D1S1()
   {
   }

   DivS1D1S1::~DivS1D1S1()
   {
   }

   void DivS1D1S1::makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const
   {
      throw std::logic_error("Operator is not implemented!");
   }

   void DivS1D1S1::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      throw std::logic_error("Operator is not implemented!");
   }

}
}
}
}
}
