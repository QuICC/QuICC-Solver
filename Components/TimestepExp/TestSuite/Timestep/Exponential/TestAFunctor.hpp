/**
 * @file TestAFunctor.hpp
 * @brief Test functor for action of matrix A
 */

#ifndef QUICC_TESTSUITE_TIMESTEP_EXPONENTIAL_TESTAFUNCTOR_HPP
#define QUICC_TESTSUITE_TIMESTEP_EXPONENTIAL_TESTAFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Timestep {

namespace Exponential {

class TestAFunctor
{
   public:
      /**
       * @brief ctor
       */
      TestAFunctor(const int n, const int id);

      /**
       * @brief ctor
       */
      TestAFunctor(const Matrix& matA);

      /**
       * @brief ctor
       */
      ~TestAFunctor() = default;

      /**
       * @brief Apply matrix A
       */
      void operator()(Eigen::Ref<Matrix> out, Eigen::Ref<Matrix> in)  const;

   private:
      /**
       * @brief Size
       */
      const int mcN;

      /**
       * @brief Matrix A
       */
      Matrix matA;
};

} // namespace Exponential
} // namespace Timestep
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_FRAMEWORK_TIMESTEP_EXPONENTIAL_TESTAFUNCTOR_HPP
