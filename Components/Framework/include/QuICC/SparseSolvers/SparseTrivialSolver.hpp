/**
 * @file SparseTrivialSolver.hpp
 * @brief Implementation of a templated (coupled) trivial solver structure
 */

#ifndef QUICC_SOLVER_SPARSETRIVIALSOLVER_HPP
#define QUICC_SOLVER_SPARSETRIVIALSOLVER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/SparseSolvers/SparseSolverBase.hpp"

namespace QuICC {

namespace Solver {

   /**
    * @brief Implementation of a templated (coupled) trivial solver structure
    */
   template <typename TOperator,typename TData,template <typename> class TSolver> class SparseTrivialSolver: public SparseSolverBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          * @param timeId  Solver timing with respect to timestepping
          * @param expTime Explicit linear timing with respect to nonlinear calculation
          */
         SparseTrivialSolver(const int start, const std::size_t timeId);

         /**
          * @brief Destructor
          */
         virtual ~SparseTrivialSolver();

         /**
          * @brief Get the number of systems in solver
          */
         std::size_t nSystem() const;

         /**
          * @brief Add RHS and solution data storage
          *
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          */
         virtual void addStorage(const int rows, const int cols);

         /**
          * @brief Initialise solution after data was copied
          */
         virtual void initSolutions();

         /**
          * @brief Update solver after updated solution was copied
          */
         virtual void updateSolutions();

         /**
          * @brief Set solver RHS data to zero
          */
         void zeroSolver();

         /**
          * @brief Set RHS data
          *
          * @param idx   Index of the data
          */
         TData& rRHSData(const int idx);

         /**
          * @brief Get solution data
          *
          * @param idx   Index of the data
          */
         const TData& solution(const int idx) const;

         /**
          * @brief Set solution data
          *
          * @param idx   Index of the data
          */
         TData& rSolution(const int idx);

         /**
          * @brief Set inhomogeneous boundary condition
          *
          * @param idx   Index of the condition
          */
         Eigen::SparseMatrix<typename TData::Scalar>& rInhomogeneous(const int idx);

      protected:
         /**
          * @brief Storage for solution of linear solve
          */
         std::vector<TData>  mSolution;

      private:
   };

   template <typename TOperator,typename TData,template <typename> class TSolver> SparseTrivialSolver<TOperator,TData,TSolver>::SparseTrivialSolver(const int start, const std::size_t timeId)
      : SparseSolverBase(start, timeId)
   {
      this->setInitialized();
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> SparseTrivialSolver<TOperator,TData,TSolver>::~SparseTrivialSolver()
   {
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseTrivialSolver<TOperator,TData,TSolver>::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add storage for solution
      this->mSolution.push_back(TData(rows,cols));
      this->mSolution.back().setZero();
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseTrivialSolver<TOperator,TData,TSolver>::initSolutions()
   {
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseTrivialSolver<TOperator,TData,TSolver>::updateSolutions()
   {
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseTrivialSolver<TOperator,TData,TSolver>::zeroSolver()
   {
      // Set solver RHS to zero
      for(unsigned int i = 0; i < this->mSolution.size(); ++i)
      {
         this->mSolution.at(i).setZero();
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> std::size_t SparseTrivialSolver<TOperator,TData,TSolver>::nSystem() const
   {
      return this->mSolution.size();
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> const TData& SparseTrivialSolver<TOperator,TData,TSolver>::solution(const int idx) const
   {
      return this->mSolution.at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> TData& SparseTrivialSolver<TOperator,TData,TSolver>::rSolution(const int idx)
   {
      return this->mSolution.at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> TData& SparseTrivialSolver<TOperator,TData,TSolver>::rRHSData(const int idx)
   {
      // WARNING: this is the same as rSolution. It's used to simplify some implementations
      return this->mSolution.at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> Eigen::SparseMatrix<typename TData::Scalar>& SparseTrivialSolver<TOperator,TData,TSolver>::rInhomogeneous(const int idx)
   {
      throw std::logic_error("Trivial solver cannot have a boundary value");
   }
}
}

#endif // QUICC_SOLVER_SPARSETRIVIALSOLVER_HPP
