/**
 * @file TestAFunctor.cpp
 * @brief Source of test functor for matrix A
 */

// System includes
//

// Project includes
//
#include "TestSuite/Timestep/Exponential/TestAFunctor.hpp"

namespace QuICC {

namespace TestSuite {

namespace Timestep {

namespace Exponential {

TestAFunctor::TestAFunctor(const int n, const int id)
   : mcN(n), matA(n, n)
{
   if(id == 0)
   {
      matA = Matrix::Random(n,n);
   }
   else
   {
      matA = Matrix::Identity(n,n);
   }
}

TestAFunctor::TestAFunctor(const Matrix& matA)
   : mcN(matA.rows()), matA(matA)
{
}

void TestAFunctor::operator()(Eigen::Ref<Matrix> out, Eigen::Ref<Matrix> in) const
{
   assert(in.rows() == this->matA.cols());
   assert(out.rows() == this->matA.rows());

   out = this->matA * in;
}

} // namespace Exponential
} // namespace Timestep
} // namespace TestSuite
} // namespace QuICC
