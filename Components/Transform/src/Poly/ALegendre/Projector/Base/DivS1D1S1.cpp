/**
 * @file DivS1D1S1.cpp
 * @brief Source of the implementation of the associated Legendre 1/Sin D Sin projector
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/DivS1D1S1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {


   void DivS1D1S1<base_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      throw std::logic_error("Operator is not implemented!");
   }

   void DivS1D1S1<base_t>::applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      throw std::logic_error("Operator is not implemented!");
   }

}
}
}
}
}
