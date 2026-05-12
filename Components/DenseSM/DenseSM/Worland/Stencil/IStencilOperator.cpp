/**
 * @file IStencilOperator.cpp
 * @brief Source of the implementation of generic dense stencil operator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "DenseSM/Worland/Stencil/IStencilOperator.hpp"
#include "Eigen/src/SparseCore/SparseUtil.h"

namespace QuICC {

namespace DenseSM {

namespace Worland {

namespace Stencil {

IStencilOperator::IStencilOperator(const int rows, const int cols,
   const Scalar_t alpha, const Scalar_t dBeta, const int l) :
    IWorlandOperator(rows, cols, alpha, dBeta), mL(l)
{
}

SparseMatrix IStencilOperator::spmat() const
{
   auto dmat = this->mat();
   SparseMatrix smat(dmat.rows(), dmat.cols());

   auto n = dmat.rows() - dmat.cols();
   assert(n > 0);
   std::vector<Eigen::Triplet<MHDFloat>> values;
   for(int j = 0; j < smat.cols(); j++)
   {
      assert(j < dmat.cols());
      for(int i = 0; i <= j + n; i++)
      {
         assert(i < dmat.rows());
         values.push_back(Eigen::Triplet<MHDFloat>(i, j, dmat(i,j)));
      }
   }

   smat.setFromTriplets(values.begin(), values.end());

   return smat;
}

} // namespace Stencil
} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
