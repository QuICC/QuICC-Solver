/**
 * @file TestAugmentedAFunctor.cpp
 * @brief Source of test functor for matrix A
 */

// System includes
//

// Project includes
//
#include "TestSuite/Timestep/Exponential/TestAugmentedAFunctor.hpp"

namespace QuICC {

namespace TestSuite {

namespace Timestep {

namespace Exponential {

TestAugmentedAFunctor::TestAugmentedAFunctor(const int n, const int id)
   : mAn(0), mBn(0), mN(n), matA(0,0), matB(0,0)
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

TestAugmentedAFunctor::TestAugmentedAFunctor(const Matrix& matA)
   : mAn(matA.rows()), mBn(0), mN(0), matA(matA), matB(0,0)
{
}

TestAugmentedAFunctor::TestAugmentedAFunctor(const Matrix& matA, const Matrix& matB)
   : mAn(matA.rows()), mBn(matB.cols()), mN(mAn + mBn), matA(matA), matB(matB)
{
}

void TestAugmentedAFunctor::operator()(Eigen::Ref<Matrix> out, Eigen::Ref<Matrix> in) const
{
   assert(in.rows() == this->mN);
   assert(out.rows() == this->mN);
   assert(out.cols() == in.cols());
   if(this->mAn == 0 || this->mBn == 0)
   {
      throw std::logic_error("Operators have not been initialized");
   }

   out.topRows(this->mAn) = this->matA * in.topRows(this->mAn) + this->matB * in.bottomRows(this->mBn);
   out.block(this->mAn, 0, this->mBn-1, out.cols()) = in.block(this->mAn + 1, 0, this->mBn-1, in.cols());
   out.bottomRows(1).array() = 0;
}

void TestAugmentedAFunctor::updateB(const Matrix& matB)
{
   this->matB = matB;
   this->mBn = matB.cols();
   this->mN = this->mAn + this->mBn;
}

} // namespace Exponential
} // namespace Timestep
} // namespace TestSuite
} // namespace QuICC
