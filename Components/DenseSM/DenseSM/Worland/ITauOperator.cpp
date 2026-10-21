/**
 * @file ITauOperator.cpp
 * @brief Source of the implementation of generic dense tau operator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "DenseSM/Worland/ITauOperator.hpp"
#include "Eigen/src/SparseCore/SparseUtil.h"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ITauOperator::ITauOperator(const int rows, const int cols,
   const Scalar_t alpha, const Scalar_t dBeta, const int l) :
    IWorlandOperator(rows, cols, alpha, dBeta), mL(l)
{
}

SparseMatrix ITauOperator::spmat() const
{
   auto dmat = this->mat();
   SparseMatrix smat(dmat.rows(), dmat.cols());

   std::vector<Eigen::Triplet<MHDFloat>> values;
   for(int j = 0; j < dmat.cols(); j++)
   {
      assert(j < dmat.cols());
      for(int i = 0; i < dmat.rows(); i++)
      {
         values.push_back(Eigen::Triplet<MHDFloat>(i, j, dmat(i,j)));
      }
   }

   smat.setFromTriplets(values.begin(), values.end());
   smat.makeCompressed();

   return smat;
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
