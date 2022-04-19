/**
 * @file SparseCoordinatorData.hpp
 * @brief Implementation of the base for a general sparse solver coordinator
 */

#ifndef QUICC_SOLVER_SPARSECOORDINATORDATA_HPP
#define QUICC_SOLVER_SPARSECOORDINATORDATA_HPP

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
#include "QuICC/TypeSelectors/TimeSchemeSelector.hpp"
#include "QuICC/Framework/Selector/ScalarField.hpp"
#include "QuICC/Framework/Selector/SparseSolver.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/SparseSolvers/SparseLinearSolver.hpp"
#include "QuICC/SparseSolvers/SparseTrivialSolver.hpp"

namespace QuICC {

namespace Solver {

   /**
    * @brief Implementation of the base for a general sparse solver coordinator
    */
   template <template <class,class,template <typename> class> class TSolver> class SparseCoordinatorData
   {
      public:
         /// Typedef for a real operator solver
         typedef TSolver<SparseMatrix,DecoupledZMatrix,Framework::Selector::SparseSolver>  RealSolverType;

         /// Typedef for a complex operator solver
         typedef TSolver<SparseMatrixZ,MatrixZ,Framework::Selector::SparseSolver>  ComplexSolverType;

         /// Typedef for a shared real operator solver
         typedef typename std::shared_ptr<RealSolverType >  SharedRealSolverType;

         /// Typedef for a shared complex operator solver
         typedef typename std::shared_ptr<ComplexSolverType >  SharedComplexSolverType;

         /// Typedef for an iterator to a real operator solver
         typedef typename std::vector<SharedRealSolverType, Eigen::aligned_allocator<SharedRealSolverType> >::iterator   RealSolver_iterator;

         /// Typedef for an iterator to a complex operator solver
         typedef typename std::vector<SharedComplexSolverType, Eigen::aligned_allocator<SharedComplexSolverType> >::iterator   ComplexSolver_iterator;

         /// Typedef for a shared scalar variable map
         typedef std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>  ScalarVariable_map;

         /// Typedef for a shared vector variable map
         typedef std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>  VectorVariable_map;

         /**
          * @brief Constructor
          */
         SparseCoordinatorData();

         /**
          * @brief Destructor
          */
         virtual ~SparseCoordinatorData();

         /**
          * @brief Set iterator for real operator solver
          *
          * \mhdBug Should not be public
          */
         void setIterator(RealSolver_iterator& solIt);

         /**
          * @brief Set iterator for complex operator solver
          *
          * \mhdBug Should not be public
          */
         void setIterator(ComplexSolver_iterator& solIt);

         /**
          * @brief Set end iterator for real operator solver
          *
          * \mhdBug Should not be public
          */
         void setEndIterator(RealSolver_iterator& solIt);

         /**
          * @brief Set end iterator for complex operator solver
          *
          * \mhdBug Should not be public
          */
         void setEndIterator(ComplexSolver_iterator& solIt);

         /**
          * @brief Get current solver time
          */
         std::size_t solveTime() const;

         /**
          * @brief Set solve time
          */
         void setSolveTime(const std::size_t timeId);

         /**
          * @brief Clear the RHS data of all solvers
          */
         void clearSolvers();

         /**
          * @brief Get current time
          */
         MHDFloat stepFraction();

      protected:
         /**
          * @brief Create a real operator
          */
         void addSolver(RealSolver_iterator solIt, const int idx, const int start, const std::size_t timeId);

         /**
          * @brief Create a complex operator
          */
         void addSolver(ComplexSolver_iterator solIt, const int idx, const int start, const std::size_t timeId);

         /**
          * @brief Initi the start rows for the solvers
          */
         void initStartRow();

         /**
          * @brief Vector of (coupled) real operator
          */
         std::vector<SharedRealSolverType, Eigen::aligned_allocator<SharedRealSolverType> > mRealSolvers;

         /**
          * @brief Vector of (coupled) complex operator
          */
         std::vector<SharedComplexSolverType, Eigen::aligned_allocator<SharedComplexSolverType> > mComplexSolvers;

         /**
          * @brief Storage for the current solve time
          */
         std::size_t   mSolveTime;

      private:
   };

   /**
    * @brief Generic implementation to setup the solver storage
    */
   template <template <class,class,template <class> class> class TSol,typename TSolverIt> void setupSolverStorage(SparseCoordinatorData<TSol>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId);

   /**
    * @brief Generic implementation to store the solver solution
    */
   template <template <class,class,template <class> class> class TSol,typename TSolverIt, typename TEq> void storeSolverSolution(SparseCoordinatorData<TSol>& coord, TEq spEq, const int idx, SpectralFieldId);

   /**
    * @brief Generic implementation to initialise the solvers
    */
   template <template <class,class,template <class> class> class TSol,typename TSolverIt> void initSolvers(SparseCoordinatorData<TSol>& coord);

   /**
    * @brief Generic implementation to update the solvers
    */
   template <template <class,class,template <class> class> class TSol,typename TSolverIt> void updateSolvers(SparseCoordinatorData<TSol>& coord);

   /**
    * @brief Generic implementation to solver solvers
    */
   template <template <class,class,template <class> class> class TSol,typename TSolverIt> std::pair<bool,MHDFloat> solveSolvers(SparseCoordinatorData<TSol>& coord);

   /**
    * @brief Generic implementation to update the time matrix solvers
    */
   template <template <class,class,template <class> class> class TSol,typename TSolverIt> void updateTimeMatrixSolvers(SparseCoordinatorData<TSol>& coord, const MHDFloat dt);

   /**
    * @brief Generic implementation to initialise the solver solution
    */
   template <template <class,class,template <class> class> class TSol,typename TSolverIt,typename TEq> void initSolverSolution(SparseCoordinatorData<TSol>& coord, TEq spEq, const int idx, SpectralFieldId id);

   /**
    * @brief Generic implementation to get the explicit linear input
    */
   template <template <class,class,template <class> class> class TSol,typename TSolverIt,typename TEq> void getExplicitSolverInput(SparseCoordinatorData<TSol>& coord, TEq spEq, const int idx, SpectralFieldId id, const std::size_t opId, const typename SparseCoordinatorData<TSol>::ScalarVariable_map& scalVar, const typename SparseCoordinatorData<TSol>::VectorVariable_map& vectVar);

   /**
    * @brief Generic implementation to get the solver input
    */
   template <template <class,class,template <class> class> class TSol,typename TSolverIt,typename TEq> void getSolverInput(SparseCoordinatorData<TSol>& coord, TEq spEq, const int idx, SpectralFieldId id, const typename SparseCoordinatorData<TSol>::ScalarVariable_map& scalVar, const typename SparseCoordinatorData<TSol>::VectorVariable_map& vectVar);

   /**
    * @brief Compute the explicit linear input independently of solver type
    */
   template <typename TSolverIt,typename TEq> void computeExplicitSolverInput(const TEq spEq, const SpectralFieldId id, const TSolverIt solveIt, const std::size_t opId, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalVar, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectVar);

   /**
    * @brief Compute the solver input independently of solver type
    */
   template <typename TSolverIt,typename TEq> void computeSolverInput(const TEq spEq, const SpectralFieldId id, const TSolverIt solveIt, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalVar, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectVar);

   /**
    * @brief Setup the storage and related information independently of solver type
    */
   template <typename TSolverIt> void setupStorage(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt);

   //
   //
   //

   template <template <class,class,template <class> class> class TSolver> SparseCoordinatorData<TSolver>::SparseCoordinatorData()
   {
   }

   template <template <class,class,template <class> class> class TSolver> SparseCoordinatorData<TSolver>::~SparseCoordinatorData()
   {
   }

   template <template <class,class,template <class> class> class TSolver> std::size_t SparseCoordinatorData<TSolver>::solveTime() const
   {
      return this->mSolveTime;
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::setSolveTime(const std::size_t timeId)
   {
      this->mSolveTime = timeId;
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::addSolver(typename SparseCoordinatorData<TSolver>::RealSolver_iterator, const int idx, const int start, const std::size_t timeId)
   {
      if(idx > static_cast<int>(this->mRealSolvers.size()) - 1)
      {
         auto spSolver = std::make_shared<SparseCoordinatorData<TSolver>::RealSolverType>(start,timeId);

         this->mRealSolvers.push_back(spSolver);
      }
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::addSolver(typename SparseCoordinatorData<TSolver>::ComplexSolver_iterator, const int idx, const int start, const std::size_t timeId)
   {
      if(idx > static_cast<int>(this->mComplexSolvers.size()) - 1)
      {
         auto spSolver = std::make_shared<SparseCoordinatorData<TSolver>::ComplexSolverType>(start,timeId);

         this->mComplexSolvers.push_back(spSolver);
      }
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::initStartRow()
   {
      for(auto rIt = this->mRealSolvers.begin(); rIt != this->mRealSolvers.end(); ++rIt)
      {
         (*rIt)->initStartRow();
      }

      for(auto zIt = this->mComplexSolvers.begin(); zIt != this->mComplexSolvers.end(); ++zIt)
      {
         (*zIt)->initStartRow();
      }
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::clearSolvers()
   {
      for(auto rIt = this->mRealSolvers.begin(); rIt != this->mRealSolvers.end(); ++rIt)
      {
         (*rIt)->zeroSolver();
      }

      for(auto zIt = this->mComplexSolvers.begin(); zIt != this->mComplexSolvers.end(); ++zIt)
      {
         (*zIt)->zeroSolver();
      }
   }

   template <template <class,class,template <class> class> class TSolver> MHDFloat SparseCoordinatorData<TSolver>::stepFraction()
   {
      if(this->mRealSolvers.begin() != this->mRealSolvers.end())
      {
         return (*this->mRealSolvers.begin())->stepFraction();
      } else
      {
         return (*this->mComplexSolvers.begin())->stepFraction();
      }
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::setIterator(typename SparseCoordinatorData<TSolver>::RealSolver_iterator& solIt)
   {
      solIt = this->mRealSolvers.begin();
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::setIterator(typename SparseCoordinatorData<TSolver>::ComplexSolver_iterator& solIt)
   {
      solIt = this->mComplexSolvers.begin();
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::setEndIterator(typename SparseCoordinatorData<TSolver>::RealSolver_iterator& solIt)
   {
      solIt = this->mRealSolvers.end();
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::setEndIterator(typename SparseCoordinatorData<TSolver>::ComplexSolver_iterator& solIt)
   {
      solIt = this->mComplexSolvers.end();
   }

   //
   //
   //
   //

   template <template <class,class,template <class> class> class TSolver, typename TSolverIt> void setupSolverStorage(SparseCoordinatorData<TSolver>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId id)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // setup storage and information
      setupStorage(spEq, id, solIt);
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt,typename TEq> void storeSolverSolution(SparseCoordinatorData<TSolver>& coord, TEq spEq, const int idx, SpectralFieldId id)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // Get solver output
      for(int i = 0; i < (*solIt)->nSystem(); i++)
      {
         spEq->storeSolution(id.second, (*solIt)->solution(i), i, (*solIt)->startRow(id,i));
      }

      // Apply constraint on solution
      spEq->applyConstraint(id.second);
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt> void initSolvers(SparseCoordinatorData<TSolver>& coord)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      TSolverIt endIt;
      coord.setEndIterator(endIt);

      for(; solIt != endIt; ++solIt)
      {
         (*solIt)->initSolver();
      }
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt> void updateSolvers(SparseCoordinatorData<TSolver>& coord)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      TSolverIt endIt;
      coord.setEndIterator(endIt);

      for(; solIt != endIt; ++solIt)
      {
         (*solIt)->updateSolver();
      }
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt> std::pair<bool,MHDFloat> solveSolvers(SparseCoordinatorData<TSolver>& coord)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      TSolverIt endIt;
      coord.setEndIterator(endIt);

      std::pair<bool,MHDFloat>  status = std::make_pair(false, -1.0);
      for(; solIt != endIt; ++solIt)
      {
         if((*solIt)->solveTiming() == coord.solveTime())
         {
            bool solving = false;
            do
            {
               // Prepare solve of linear system
               bool needSolve = (*solIt)->preSolve();

               if(needSolve)
               {
                  // Solve linear system
                  (*solIt)->solve();

                  // Work on fields after solve
                  solving = (*solIt)->postSolve();

               } else
               {
                  solving = false;
               }

            } while (solving);

            status.first = (*solIt)->finished();

            if(status.first)
            {
               status.second = std::max(status.second, (*solIt)->error());
            }
         }
      }

      return status;
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt> void updateTimeMatrixSolvers(SparseCoordinatorData<TSolver>& coord, const MHDFloat dt)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      TSolverIt endIt;
      coord.setEndIterator(endIt);

      for(; solIt != endIt; ++solIt)
      {
         // Compute linear solve RHS
         (*solIt)->updateTimeMatrix(dt);
      }
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt,typename TEq> void initSolverSolution(SparseCoordinatorData<TSolver>& coord, TEq spEq, const int idx, SpectralFieldId id)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // Get solver input
      for(int i = 0; i < (*solIt)->nSystem(); i++)
      {
         if(spEq->couplingInfo(id.second).isGalerkin())
         {
            Equations::solveStencilUnknown(*spEq, id.second, (*solIt)->rSolution(i), i, (*solIt)->startRow(id,i));

         } else
         {
            std::visit([&](auto&& p){Equations::copyUnknown(*spEq, p->dom(0).perturbation(), id.second, (*solIt)->rSolution(i), i, (*solIt)->startRow(id,i), true, true);}, spEq->spUnknown());
         }
      }

      (*solIt)->initSolutions();
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt,typename TEq> void getExplicitSolverInput(SparseCoordinatorData<TSolver>& coord, TEq spEq, const int idx, SpectralFieldId id, const std::size_t opId, const typename SparseCoordinatorData<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorData<TSolver>::VectorVariable_map& vectVar)
   {
      // Create iterator to current complex field solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // Get solver input
      computeExplicitSolverInput<TSolverIt,TEq>(spEq, id, solIt, opId, scalVar, vectVar);
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt,typename TEq> void getSolverInput(SparseCoordinatorData<TSolver>& coord, TEq spEq, const int idx, SpectralFieldId id, const typename SparseCoordinatorData<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorData<TSolver>::VectorVariable_map& vectVar)
   {
      // Create iterator to current complex field solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // Get solver input
      computeSolverInput<TSolverIt,TEq>(spEq, id, solIt, scalVar, vectVar);
   }

   //
   //
   //

   template <typename TSolverIt,typename TEq> void computeExplicitSolverInput(const TEq spEq, const SpectralFieldId id, const TSolverIt solveIt, const std::size_t opId, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalVar, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectVar)
   {
      // Get timestep input
      for(int i = 0; i < (*solveIt)->nSystem(); i++)
      {
         // Loop over explicit fields
         for(auto& fIt: make_range(spEq->couplingInfo(id.second).explicitRange(opId)))
         {
            if(fIt.second == FieldComponents::Spectral::SCALAR)
            {
               std::visit([&](auto&& p){Equations::addExplicitTerm(*spEq, opId, id.second, (*solveIt)->rRHSData(i), (*solveIt)->startRow(id,i), fIt, p->dom(0).perturbation(), i);}, scalVar.find(fIt.first)->second);
            } else
            {
               std::visit([&](auto&& p){Equations::addExplicitTerm(*spEq, opId, id.second, (*solveIt)->rRHSData(i), (*solveIt)->startRow(id,i), fIt, p->dom(0).perturbation().comp(fIt.second), i);}, vectVar.find(fIt.first)->second);
            }
         }
      }
   }

   template <typename TSolverIt,typename TEq> void computeSolverInput(const TEq spEq, const SpectralFieldId id, const TSolverIt solveIt, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>&, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>&)
   {
      // Get timestep input
      for(int i = 0; i < (*solveIt)->nSystem(); i++)
      {
         // Copy field values into solver input
         Equations::copyNonlinear(*spEq, id.second, (*solveIt)->rRHSData(i), i, (*solveIt)->startRow(id,i));

         // Add source term
         std::visit([&](auto&& p){Equations::addSource(*spEq, p->dom(0).perturbation(), id.second, (*solveIt)->rRHSData(i), i, (*solveIt)->startRow(id,i));}, spEq->spUnknown());

         // If required set inhomogenous boundary condition value
         if(spEq->couplingInfo(id.second).hasBoundaryValue())
         {
            // Set boundary value
            std::visit([&](auto&& p){Equations::setBoundaryValue(*spEq, p->dom(0).perturbation(), id.second, (*solveIt)->rInhomogeneous(i), i, (*solveIt)->startRow(id,i));}, spEq->spUnknown());
         }
      }
   }

   template <typename TSolverIt> void setupStorage(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt)
   {
      // Number of linear systems
      int nSystems = spEq->couplingInfo(id.second).nSystems();

      // Start row for storage information
      ArrayI startRow(nSystems);

      // Initialise the linear solver
      if((*solveIt)->nSystem() == 0)
      {
         // Initialise field storage and information
         for(int i = 0; i < nSystems; i++)
         {
            // Create data storage
            (*solveIt)->addStorage(spEq->couplingInfo(id.second).systemN(i), spEq->couplingInfo(id.second).rhsCols(i));
         }
      }

      // Gather information about the solver matrices
      for(int i = 0; i < nSystems; i++)
      {
         // Store the start row
         startRow(i) = spEq->couplingInfo(id.second).galerkinN(i);
      }

      // Store storage information
      (*solveIt)->addInformation(id, spEq->couplingInfo(id.second).fieldIndex(), startRow);
   }

}
}

#endif // QUICC_SOLVER_SPARSECOORDINATORDATA_HPP
