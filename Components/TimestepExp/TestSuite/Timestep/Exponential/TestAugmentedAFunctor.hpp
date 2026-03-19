/**
 * @file TestAugmentedAFunctor.hpp
 * @brief Test functor for action of matrix A
 */

#ifndef QUICC_TESTSUITE_TIMESTEP_EXPONENTIAL_TESTAUGMENTEDFUNCTOR_HPP
#define QUICC_TESTSUITE_TIMESTEP_EXPONENTIAL_TESTAUGMENTEDFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Timestep {

namespace Exponential {

class TestAugmentedAFunctor
{
   public:
      /**
       * @brief ctor
       */
      TestAugmentedAFunctor(const int n, const int id);

      /**
       * @brief ctor
       */
      TestAugmentedAFunctor(const Matrix& matA);

      /**
       * @brief ctor
       */
      TestAugmentedAFunctor(const Matrix& matA, const Matrix& matB);

      /**
       * @brief ctor
       */
      ~TestAugmentedAFunctor() = default;

      /**
       * @brief Apply matrix A
       */
      void operator()(Eigen::Ref<Matrix> out, Eigen::Ref<Matrix> in)  const;

      /**
       * @brief Update matrix B
       */
      void updateB(const Matrix &matB);

   private:
      /**
       * @brief Size if A
       */
      int mAn;

      /**
       * @brief Size if B
       */
      int mBn;

      /**
       * @brief Total size
       */
      int mN;

      /**
       * @brief Matrix A
       */
      Matrix matA;

      /**
       * @brief Matrix B
       */
      Matrix matB;
};

} // namespace Exponential
} // namespace Timestep
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_TIMESTEP_EXPONENTIAL_TESTAUGMENTEDFUNCTOR_HPP
