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

TestAFunctor::TestAFunctor()
{}

void TestAFunctor::operator()(Eigen::Ref<Matrix> out, Eigen::Ref<Matrix> in) const
{
   Matrix matA = Matrix::Random(in.rows(), in.cols());

   out = matA * in;
}

} // namespace Exponential
} // namespace Timestep
} // namespace TestSuite
} // namespace QuICC
