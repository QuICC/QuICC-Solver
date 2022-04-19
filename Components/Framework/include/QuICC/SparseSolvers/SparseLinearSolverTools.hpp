/**
 * @file SparseLinearSolverTools.hpp
 * @brief Implementation of a couple of templated linear solver helper functions
 */

#ifndef QUICC_SOLVER_SPARSELINEARSOLVERTOOLS_HPP
#define QUICC_SOLVER_SPARSELINEARSOLVERTOOLS_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace Solver {

namespace internal {

   /**
    * @brief Simpler wrapper around solver solve call
    */
   template <typename TSolver> void solveWrapper(Eigen::Matrix<typename TSolver::MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>& rSolution, TSolver& solver, const Eigen::Matrix<typename TSolver::MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>& rhs);
   template <typename TSolver> void solveWrapper(Eigen::Matrix<std::complex<typename TSolver::MatrixType::Scalar>, Eigen::Dynamic, Eigen::Dynamic>& rSolution, TSolver& solver, const Eigen::Matrix<std::complex<typename TSolver::MatrixType::Scalar>, Eigen::Dynamic, Eigen::Dynamic>& rhs);
   template <typename TSolver> void solveWrapper(DecoupledZMatrix& rSolution, TSolver& solver, const DecoupledZMatrix& rhs);

   /**
    * @brief Simpler wrapper around solver solve call for shared pointer to solver
    */
   template <typename TSolver> void solveWrapper(Eigen::Matrix<typename TSolver::MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>& rSolution, std::shared_ptr<TSolver> spSolver, const Eigen::Matrix<typename TSolver::MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>& rhs);
   template <typename TSolver> void solveWrapper(Eigen::Matrix<std::complex<typename TSolver::MatrixType::Scalar>, Eigen::Dynamic, Eigen::Dynamic>& rSolution, std::shared_ptr<TSolver> spSolver, const Eigen::Matrix<std::complex<typename TSolver::MatrixType::Scalar>, Eigen::Dynamic, Eigen::Dynamic>& rhs);
   template <typename TSolver> void solveWrapper(DecoupledZMatrix& rSolution, std::shared_ptr<TSolver> spSolver, const DecoupledZMatrix& rhs);
}

namespace internal {

   template <typename TSolver>  inline void solveWrapper(Eigen::Matrix<typename TSolver::MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>& rSolution, TSolver& solver, const Eigen::Matrix<typename TSolver::MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>& rhs)
   {
      rSolution = solver.solve(rhs);

      // Safety assert for successful solve
      assert(solver.info() == Eigen::Success);
   }

   template <typename TSolver> inline void solveWrapper(Eigen::Matrix<std::complex<typename TSolver::MatrixType::Scalar>, Eigen::Dynamic, Eigen::Dynamic>& rSolution, TSolver& solver, const Eigen::Matrix<std::complex<typename TSolver::MatrixType::Scalar>, Eigen::Dynamic, Eigen::Dynamic>& rhs)
   {
      Matrix tmpIn(rhs.rows(), rhs.cols());
      Matrix tmpOut(rhs.rows(), rhs.cols());

      tmpIn = rhs.real();
      tmpOut = solver.solve(tmpIn);
      rSolution.real() = tmpOut;

      // Safety assert for successful solve
      assert(solver.info() == Eigen::Success);

      tmpIn = rhs.imag();
      tmpOut = solver.solve(tmpIn);
      rSolution.imag() = tmpOut;

      // Safety assert for successful solve
      assert(solver.info() == Eigen::Success);
   }

   template <typename TSolver> inline void solveWrapper(DecoupledZMatrix& rSolution, TSolver& solver, const DecoupledZMatrix& rhs)
   {
      rSolution.real() = solver.solve(rhs.real());

      // Safety assert for successful solve
      assert(solver.info() == Eigen::Success);

      rSolution.imag() = solver.solve(rhs.imag());

      // Safety assert for successful solve
      assert(solver.info() == Eigen::Success);
   }

   template <typename TSolver> inline void solveWrapper(Eigen::Matrix<typename TSolver::MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>& rSolution, std::shared_ptr<TSolver> spSolver, const Eigen::Matrix<typename TSolver::MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>& rhs)
   {
      rSolution = spSolver->solve(rhs);

      // Safety assert for successful solve
      assert(spSolver->info() == Eigen::Success);
   }

   template <typename TSolver> inline void solveWrapper(Eigen::Matrix<std::complex<typename TSolver::MatrixType::Scalar>, Eigen::Dynamic, Eigen::Dynamic>& rSolution, std::shared_ptr<TSolver> spSolver, const Eigen::Matrix<std::complex<typename TSolver::MatrixType::Scalar>, Eigen::Dynamic, Eigen::Dynamic>& rhs)
   {
      Matrix tmp(rhs.rows(), rhs.cols());
      tmp = rhs.real();
      rSolution.real() = spSolver->solve(tmp);

      // Safety assert for successful solve
      assert(spSolver->info() == Eigen::Success);

      tmp = rhs.imag();
      rSolution.imag() = spSolver->solve(tmp);

      // Safety assert for successful solve
      assert(spSolver->info() == Eigen::Success);
   }

   template <typename TSolver> inline void solveWrapper(DecoupledZMatrix& rSolution, std::shared_ptr<TSolver> spSolver, const DecoupledZMatrix& rhs)
   {
      rSolution.real() = spSolver->solve(rhs.real());

      // Safety assert for successful solve
      assert(spSolver->info() == Eigen::Success);

      rSolution.imag() = spSolver->solve(rhs.imag());

      // Safety assert for successful solve
      assert(spSolver->info() == Eigen::Success);
   }
}
}
}

#endif // QUICC_SOLVER_SPARSELINEARSOLVERTOOLS_HPP
