/**
 * @file TestHelper.hpp
 * @brief Helper functions to setup Timestep tests
 */

#ifndef QUICC_TESTSUITE_TIMESTEP_RUNGEKUTTACB_TESTHELPER_HPP
#define QUICC_TESTSUITE_TIMESTEP_RUNGEKUTTACB_TESTHELPER_HPP

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/Solver/SparseSolver.hpp"
#include "TestSuite/Timestep/RungeKuttaCB/Test.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Timestep {

namespace RungeKuttaCB {

void createOperators(Test& test, std::map<std::size_t, DecoupledZSparse>& ops);
void initQuasiInverse(Test& test, SparseMatrix& qi);
void setForcing(Test& test, Matrix& forcing, const MHDFloat t);
void initSolution(Test& test, DecoupledZMatrix& sol);
void checkSolution(Test& test, const Matrix& data, const int tIdx);

template <typename TScheme>
std::shared_ptr<typename TScheme::template ImplementationType<SparseMatrix,
   DecoupledZMatrix, ::QuICC::Framework::Selector::SparseSolver>>
createTimestepper(Test& test)
{
   typedef typename TScheme::template ImplementationType<SparseMatrix,
      DecoupledZMatrix, ::QuICC::Framework::Selector::SparseSolver>
      Timestepper;

   const auto& start = test.start;
   auto spStepper =
      std::make_shared<Timestepper>(start, SolveTiming::Prognostic::id());

   auto spScheme = std::make_shared<TScheme>();
   spStepper->setScheme(spScheme);

   return spStepper;
}

template <typename TScheme>
void setupTimestepper(Test& test,
   std::shared_ptr<typename TScheme::template ImplementationType<SparseMatrix,
      DecoupledZMatrix, ::QuICC::Framework::Selector::SparseSolver>>
      spStepper)
{
   const int nSystems = 1;
   const auto& nN = test.nN;
   const auto& cols = test.cols;
   const auto& startRow = test.startRow;
   const auto& dt = test.dt;
   const int idx = 0;

   spStepper->addStorage(nN, cols);
   auto fieldId = std::make_pair(PhysicalNames::Temperature::id(),
      FieldComponents::Spectral::SCALAR);
   spStepper->addInformation(fieldId, 0, startRow);
   spStepper->initStartRow();

   spStepper->initMatrices(nSystems);

   // Build operators
   std::map<std::size_t, DecoupledZSparse> ops;
   createOperators(test, ops);
   spStepper->buildOperators(idx, ops, dt, nN);

   // Solver is initialized
   spStepper->setInitialized();
   spStepper->initSolver();

   // Create initial state
   initSolution(test, spStepper->rSolution(idx));
   spStepper->initSolutions();
}

template <typename TScheme>
void singleStepForward(Test& test,
   std::shared_ptr<typename TScheme::template ImplementationType<SparseMatrix,
      DecoupledZMatrix, ::QuICC::Framework::Selector::SparseSolver>>
      spStepper)
{
   const auto& nN = test.nN;
   const auto& cols = test.cols;
   const auto& tEnd = test.tEnd;
   const int idx = 0;

   // Initialize quasi-inverse operators
   SparseMatrix qi;
   initQuasiInverse(test, qi);

   std::pair<bool, MHDFloat> status = std::make_pair(false, -1.0);

   // Test single steps
   for (int tIdx = 0; tIdx < tEnd.size(); tIdx++)
   {
      MHDFloat t = test.t0;
      auto t_ = t;
      // Update timestep
      const auto& dt = tEnd(tIdx);
      spStepper->updateTimeMatrix(dt);
      spStepper->updateSolver();

      // Reset initial state
      initSolution(test, spStepper->rSolution(idx));
      spStepper->initSolutions();

      do
      {
         // Set RHS
         Matrix forcing(nN, cols);
         setForcing(test, forcing, t_);
         spStepper->rRHSData(idx).real() = qi * forcing;
         spStepper->rRHSData(idx).imag() = qi * forcing;

         bool solving = false;
         do
         {
            // Prepare solve of linear system
            bool needSolve = spStepper->preSolve();

            if (needSolve)
            {
               // Solve linear system
               spStepper->solve();

               // Work on fields after solve
               solving = spStepper->postSolve();
            }
            else
            {
               solving = false;
            }

         } while (solving);

         status.first = spStepper->finished();

         if (status.first)
         {
            status.second = std::max(status.second, spStepper->error());
            t += dt;
            t_ = t;
         }
         else
         {
            t_ = t + spStepper->stepFraction() * dt;
         }
      } while (!status.first);

      // Check solution
      checkSolution(test, spStepper->solution(idx).real(), tIdx);
   }
}

template <typename TScheme>
void multiStepForward(Test& test,
   std::shared_ptr<typename TScheme::template ImplementationType<SparseMatrix,
      DecoupledZMatrix, ::QuICC::Framework::Selector::SparseSolver>>
      spStepper,
   const int steps)
{
   const auto& nN = test.nN;
   const auto& cols = test.cols;
   const auto& tEnd = test.tEnd;
   const int idx = 0;

   // Initialize quasi-inverse operators
   SparseMatrix qi;
   initQuasiInverse(test, qi);

   std::pair<bool, MHDFloat> status = std::make_pair(false, -1.0);

   // Test single steps
   for (int tIdx = 0; tIdx < tEnd.size(); tIdx++)
   {
      MHDFloat t = test.t0;
      auto t_ = t;
      // Update timestep
      const auto& dt = tEnd(tIdx) / static_cast<MHDFloat>(steps);
      spStepper->updateTimeMatrix(dt);
      spStepper->updateSolver();

      // Reset initial state
      initSolution(test, spStepper->rSolution(idx));
      spStepper->initSolutions();

      do
      {
         do
         {
            // Set RHS
            Matrix forcing(nN, cols);
            setForcing(test, forcing, t_);
            spStepper->rRHSData(idx).real() = qi * forcing;
            spStepper->rRHSData(idx).imag() = qi * forcing;

            bool solving = false;
            do
            {
               // Prepare solve of linear system
               bool needSolve = spStepper->preSolve();

               if (needSolve)
               {
                  // Solve linear system
                  spStepper->solve();

                  // Work on fields after solve
                  solving = spStepper->postSolve();
               }
               else
               {
                  solving = false;
               }

            } while (solving);

            status.first = spStepper->finished();

            if (status.first)
            {
               status.second = std::max(status.second, spStepper->error());
               t += dt;
               t_ = t;
            }
            else
            {
               t_ = t + spStepper->stepFraction() * dt;
            }
         } while (!status.first);
      } while (t < tEnd(tIdx));

      // Check solution
      checkSolution(test, spStepper->solution(idx).real(), tIdx);
   }
}

} // namespace RungeKuttaCB
} // namespace Timestep
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_TIMESTEP_RUNGEKUTTACB_TESTHELPER_HPP
