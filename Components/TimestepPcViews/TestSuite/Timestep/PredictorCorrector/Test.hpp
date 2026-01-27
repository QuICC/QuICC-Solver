/**
 * @file TestHelper.hpp
 * @brief Helper functions to setup Timestep tests
 */

#ifndef QUICC_TESTSUITE_TIMESTEP_PREDICTORCORRECTOR_TEST_HPP
#define QUICC_TESTSUITE_TIMESTEP_PREDICTORCORRECTOR_TEST_HPP

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/Solver/SparseSolver.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Timestep {

namespace PredictorCorrector {

class Test
{
public:
   enum class BasisId
   {
      WORLAND,
      CHEBYSHEV
   };

   enum class EquationId
   {
      DIFFUSION,
      BIDIFFUSION
   };

   Test(const int id);
   ~Test() = default;

   /**
    * @brief Get reference data
    */
   Matrix getReference();

   /**
    * @brief Get initial state
    */
   Matrix getInitial();

   /**
    * @brief Get forcing data
    */
   Matrix getForcing();

   int nN;

   /**
    * @brief Basis ID
    */
   BasisId basisId;

   /**
    * @brief Equation ID
    */
   EquationId equationId;

   int l;
   MHDFloat dt;
   Array timeCoeffs;
   Array tEnd;
   ArrayI startRow;

   /**
    * @brief Filename base
    */
   std::string fbase;

   /**
    * @brief Data start index
    */
   const int start;

   /**
    * @brief number of colums
    */
   const int cols;

   /**
    * @brief Starting time
    */
   const MHDFloat t0;

private:
   /**
    * @brief Configure test
    */
   void configure();

   /**
    * @brief ID of tests
    */
   int id;
};

} // namespace PredictorCorrector
} // namespace Timestep
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_TIMESTEP_PREDICTORCORRECTOR_TEST_HPP
