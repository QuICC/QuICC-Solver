/**
 * @file SparseDummySolver.hpp
 * @brief Implementation of the sparse dummy solver
 */

#ifndef QUICC_SOLVER_SPARSEDUMMYSOLVER_HPP
#define QUICC_SOLVER_SPARSEDUMMYSOLVER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Solver/SparseSolver.hpp"

namespace QuICC {

namespace Solver {

   /**
    * @brief Implementation of the sparse dummy solver
    */
   template <typename TOperator, typename TData, template <typename> class TSolver> class SparseDummySolver
   {
      public:
         /**
          * @brief Dummy constructor
          *
          * @param start   Starting index (for example without m=0)
          * @param time    Solver timing with respect to timestepping
          */
         SparseDummySolver(const int start, const std::size_t timeId) {};

         /**
          * @brief Dummy destructor
          */
         virtual ~SparseDummySolver() {};

         /**
          * @brief Dummy implementation
          */
         void initStartRow();

         /**
          * @brief Dummy implementation
          */
         void zeroSolver();

         /**
          * @brief Get current timestep fraction
          */
         MHDFloat stepFraction() const;

      protected:

      private:
   };

   template <typename TOperator, typename TData, template <typename> class TSolver> void SparseDummySolver<TOperator,TData,TSolver>::initStartRow()
   {
      throw std::logic_error("Dummy solver should not be used!");
   }

   template <typename TOperator, typename TData, template <typename> class TSolver> void SparseDummySolver<TOperator,TData,TSolver>::zeroSolver()
   {
      throw std::logic_error("Dummy solver should not be used!");
   }

   template <typename TOperator, typename TData, template <typename> class TSolver> MHDFloat SparseDummySolver<TOperator,TData,TSolver>::stepFraction() const
   {
      throw std::logic_error("Dummy solver should not be used!");
   }

   typedef SparseDummySolver<SparseMatrixZ,MatrixZ,Framework::Selector::SparseSolver> SparseDummySolverComplexType;

   typedef std::shared_ptr<SparseDummySolverComplexType> SharedSparseDummySolverComplexType;

   typedef std::vector<SharedSparseDummySolverComplexType>::iterator ComplexDummy_iterator;

} // Solver
} // QuICC

#endif // QUICC_SOLVER_SPARSEDUMMYSOLVER_HPP
