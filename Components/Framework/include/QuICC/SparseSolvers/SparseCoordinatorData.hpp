/**
 * @file SparseCoordinatorData.hpp
 * @brief Implementation of the base for a general sparse solver coordinator
 */

#ifndef QUICC_SOLVER_SPARSECOORDINATORDATA_HPP
#define QUICC_SOLVER_SPARSECOORDINATORDATA_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "Environment/QuICCEnv.hpp"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/ModelOperator/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#endif
#include "Types/Typedefs.hpp"
#include "QuICC/SolveTiming/Before.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Solver/SparseSolver.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/Equations/StoreSolution.hpp"
#include "QuICC/Equations/AddSource.hpp"
#include "QuICC/Equations/SetBoundaryValue.hpp"
#include "QuICC/Equations/ExplicitTerm.hpp"
#include "QuICC/Equations/CopyNonlinear.hpp"
#include "QuICC/Equations/SolveStencilUnknown.hpp"
#include "QuICC/Equations/CorrectSolution.hpp"
#include "QuICC/SparseSolvers/SparseLinearSolver.hpp"
#include "QuICC/SparseSolvers/SparseTrivialSolver.hpp"
#include "ViewOps/Transpose/OpGrouped.hpp"

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
         typedef typename std::vector<SharedRealSolverType>::iterator   RealSolver_iterator;

         /// Typedef for an iterator to a complex operator solver
         typedef typename std::vector<SharedComplexSolverType>::iterator   ComplexSolver_iterator;

         /// Typedef for a shared scalar variable map
         typedef std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>  ScalarVariable_map;

         /// Typedef for a shared vector variable map
         typedef std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>  VectorVariable_map;

         /**
          * @brief Constructor
          */
         SparseCoordinatorData() = default;

         /**
          * @brief Destructor
          */
         virtual ~SparseCoordinatorData() = default;

         /**
          * @brief Set iterator for real operator solver
          *
          * @param solIt   Solver iterator
          * @param idx     Index
          */
         void setIterator(RealSolver_iterator& solIt, const int idx = 0);

         /**
          * @brief Set iterator for complex operator solver
          *
          * @param solIt   Solver iterator
          * @param idx     Index
          */
         void setIterator(ComplexSolver_iterator& solIt, const int idx = 0);

         /**
          * @brief Set iterator range for real operator solver
          *
          * @param begIt   Begin iterator
          * @param endIt   End iterator
          */
         void setRange(RealSolver_iterator& begIt, RealSolver_iterator& endIt);

         /**
          * @brief Set iterator range for complex operator solver
          *
          * @param begIt   Begin iterator
          * @param endIt   End iterator
          */
         void setRange(ComplexSolver_iterator& begIt, ComplexSolver_iterator& endIt);

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
         std::vector<SharedRealSolverType> mRealSolvers;

         /**
          * @brief Vector of (coupled) complex operator
          */
         std::vector<SharedComplexSolverType> mComplexSolvers;

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
      DebuggerMacro_msg("Creating real solver for idx = " + std::to_string(idx) + "/" + std::to_string(this->mRealSolvers.size()), 3);

      for(int i = static_cast<int>(this->mRealSolvers.size()); i <= idx; i++)
      {
         this->mRealSolvers.push_back(nullptr);
      }

      if(this->mRealSolvers.at(idx) == nullptr)
      {
         auto spSolver = std::make_shared<SparseCoordinatorData<TSolver>::RealSolverType>(start,timeId);

         this->mRealSolvers.at(idx) = spSolver;
      }

      DebuggerMacro_msg("... current number of real solvers = " + std::to_string(this->mRealSolvers.size()), 3);
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::addSolver(typename SparseCoordinatorData<TSolver>::ComplexSolver_iterator, const int idx, const int start, const std::size_t timeId)
   {
      DebuggerMacro_msg("Creating complex solver for idx = " + std::to_string(idx) + "/" + std::to_string(this->mComplexSolvers.size()), 3);

      for(int i = static_cast<int>(this->mComplexSolvers.size()); i <= idx; i++)
      {
         this->mComplexSolvers.push_back(nullptr);
      }

      if(this->mComplexSolvers.at(idx) == nullptr)
      {
         auto spSolver = std::make_shared<SparseCoordinatorData<TSolver>::ComplexSolverType>(start,timeId);

         this->mComplexSolvers.at(idx) = spSolver;
      }

      DebuggerMacro_msg("... current number of complex solvers = " + std::to_string(this->mComplexSolvers.size()), 3);
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
         if((*rIt)->solveTiming() == this->solveTime())
         {
            (*rIt)->zeroSolver();
         }
      }

      for(auto zIt = this->mComplexSolvers.begin(); zIt != this->mComplexSolvers.end(); ++zIt)
      {
         if((*zIt)->solveTiming() == this->solveTime())
         {
            (*zIt)->zeroSolver();
         }
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

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::setIterator(typename SparseCoordinatorData<TSolver>::RealSolver_iterator& solIt, const int idx)
   {
      solIt = this->mRealSolvers.begin();
      assert(std::distance(solIt, this->mRealSolvers.end()) > idx);
      std::advance(solIt, idx);
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::setIterator(typename SparseCoordinatorData<TSolver>::ComplexSolver_iterator& solIt, const int idx)
   {
      solIt = this->mComplexSolvers.begin();
      assert(std::distance(solIt, this->mComplexSolvers.end()) > idx);
      std::advance(solIt, idx);
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::setRange(typename SparseCoordinatorData<TSolver>::RealSolver_iterator& begIt, typename SparseCoordinatorData<TSolver>::RealSolver_iterator& endIt)
   {
      begIt = this->mRealSolvers.begin();
      endIt = this->mRealSolvers.end();
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorData<TSolver>::setRange(typename SparseCoordinatorData<TSolver>::ComplexSolver_iterator& begIt, typename SparseCoordinatorData<TSolver>::ComplexSolver_iterator& endIt)
   {
      begIt = this->mComplexSolvers.begin();
      endIt = this->mComplexSolvers.end();
   }

   //
   //
   //
   //

   template <template <class,class,template <class> class> class TSolver, typename TSolverIt> void setupSolverStorage(SparseCoordinatorData<TSolver>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId id)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt, idx);

      // setup storage and information
      setupStorage(spEq, id, solIt);
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt,typename TEq> void storeSolverSolution(SparseCoordinatorData<TSolver>& coord, TEq spEq, const int idx, SpectralFieldId id)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt, idx);

      DebuggerMacro_msg("Get solver solution for " + PhysicalNames::Coordinator::tag(id.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(id.second)) + ")", 6);

      // Get solver output
#ifdef QUICC_DEBUG
      spEq->corruptUnknown(id.second);
#endif //QUICC_DEBUG
      const auto& cinfo = spEq->couplingInfo(id.second);
      const SparseMatrix* pOp;
      for(std::size_t i = 0; i < (*solIt)->nSystem(); i++)
      {
         if(cinfo.isGalerkin())
         {
            pOp = &spEq->galerkinStencils(id.second).at(i);
         }
         std::visit(
               [&](auto&& p)
               {
                  Equations::storeSolution(p->rDom(0).rPerturbation().rComp(id.second), spEq->res(), cinfo, pOp, spEq->spSolutionUpdater(id.second), (*solIt)->solution(i), i, (*solIt)->startRow(id,i));
               }, spEq->spUnknown());
      }

      // Apply constraint on solution
      auto changedSolution = spEq->applyConstraint(id.second, SolveTiming::After::id());

      // Update solver solution if constraint modified it
      if(changedSolution)
      {
         auto corr = spEq->correctionConstraint(id.second, SolveTiming::After::id());
         for(std::size_t i = 0; i < (*solIt)->nSystem(); i++)
         {
            Equations::correctSolution(spEq->res(), cinfo, corr, (*solIt)->rSolution(i), i, (*solIt)->startRow(id,i));
         }

         (*solIt)->updateSolutions();
      }
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt> void initSolvers(SparseCoordinatorData<TSolver>& coord)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      TSolverIt endIt;
      coord.setRange(solIt, endIt);

      for(; solIt != endIt; ++solIt)
      {
         (*solIt)->initSolver();
      }
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt> void updateSolvers(SparseCoordinatorData<TSolver>& coord)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      TSolverIt endIt;
      coord.setRange(solIt, endIt);

      for(; solIt != endIt; ++solIt)
      {
         (*solIt)->updateSolver();
      }
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt> std::pair<bool,MHDFloat> solveSolvers(SparseCoordinatorData<TSolver>& coord)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      TSolverIt endIt;
      coord.setRange(solIt, endIt);

      std::pair<bool,MHDFloat>  status = std::make_pair(false, -1.0);
      for(; solIt != endIt; ++solIt)
      {
         if((*solIt)->solveTiming() == coord.solveTime())
         {
            bool solving = false;
            DebuggerMacro_msg("Looping over Solvers", 8);
            do
            {
               // Prepare solve of linear system
               bool needSolve = (*solIt)->preSolve();

               if(needSolve)
               {
                  DebuggerMacro_msg("...solving", 9);

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
      TSolverIt endIt;
      coord.setRange(solIt, endIt);

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
      coord.setIterator(solIt, idx);

      // Get solver input
      for(std::size_t i = 0; i < (*solIt)->nSystem(); i++)
      {
         if(spEq->couplingInfo(id.second).isGalerkin())
         {
            std::visit(
                  [&](auto&& p)
                  {
                     Equations::solveStencilUnknown(spEq->res(), spEq->couplingInfo(id.second), id, p->dom(0).perturbation().comp(id.second), (*solIt)->rSolution(i), i, (*solIt)->startRow(id,i), spEq->backend(), spEq->bcIds().map(), spEq->eqParams().map());
                  }, spEq->spUnknown());
         }
         else
         {
            std::visit(
                  [&](auto&& p)
                  {
                     Equations::copyUnknown(spEq->res(), spEq->couplingInfo(id.second), p->dom(0).perturbation().comp(id.second), (*solIt)->rSolution(i), i, (*solIt)->startRow(id,i), true, true, true);
                  }, spEq->spUnknown());
         }
      }

      (*solIt)->initSolutions();
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt,typename TEq> void getExplicitSolverInput(SparseCoordinatorData<TSolver>& coord, TEq spEq, const int idx, SpectralFieldId id, const std::size_t opId, const typename SparseCoordinatorData<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorData<TSolver>::VectorVariable_map& vectVar)
   {
      // Create iterator to current complex field solver
      TSolverIt solIt;
      coord.setIterator(solIt, idx);

      DebuggerMacro_msg("Get explicit solver input for " + PhysicalNames::Coordinator::tag(id.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(id.second)) + ")", 6);

      // Get solver input
      computeExplicitSolverInput<TSolverIt,TEq>(spEq, id, solIt, opId, scalVar, vectVar);
   }

   template <template <class,class,template <class> class> class TSolver,typename TSolverIt,typename TEq> void getSolverInput(SparseCoordinatorData<TSolver>& coord, TEq spEq, const int idx, SpectralFieldId id, const typename SparseCoordinatorData<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorData<TSolver>::VectorVariable_map& vectVar)
   {
      // Apply constraint on solution
      spEq->applyConstraint(id.second, SolveTiming::Before::id());

      // Create iterator to current complex field solver
      TSolverIt solIt;
      coord.setIterator(solIt, idx);

      DebuggerMacro_msg("Get solver input for " + PhysicalNames::Coordinator::tag(id.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(id.second)) + ")", 6);

      // Get solver input
      computeSolverInput<TSolverIt,TEq>(spEq, id, solIt, scalVar, vectVar);
   }

   //
   //
   //

   template <typename TSolverIt,typename TEq> void computeExplicitSolverInput(const TEq spEq, const SpectralFieldId id, const TSolverIt solveIt, const std::size_t opId, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalVar, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectVar)
   {
      // Build range of operator
      auto r = make_range(spEq->couplingInfo(id.second).explicitRange(opId));

#ifdef QUICC_DEBUG
      if(r.size() == 0)
      {
         DebuggerMacro_msg("(Nothing)", 7);
      }
#endif // QUICC_DEBUG

      auto addExp = [&](auto&& spEq, auto&& fIt, auto&& i, auto&& field)
      {
         std::variant<const SparseMatrix*,const SparseMatrixZ*> pOp;

         // Compute with complex linear operator
         if (spEq->hasExplicitZTerm(opId, id.second, fIt))
         {
            // Create pointer to sparse operator
            pOp = &spEq->template explicitOperators<SparseMatrixZ>(opId,
                  id.second, fIt).at(i);
         }

         // Compute with real linear operator
         if (spEq->hasExplicitDTerm(opId, id.second, fIt))
         {
            // Create pointer to sparse operator
            pOp = &spEq->template explicitOperators<SparseMatrix>(opId,
                  id.second, fIt).at(i);
         }

         std::visit(
               [&](auto& p)
               {
                  Equations::addExplicitTerm(spEq->res(), spEq->couplingInfo(id.second), *p, (*solveIt)->rRHSData(i), (*solveIt)->startRow(id,i), field, i, true, false);
               }, pOp);
      };

      // Loop over explicit fields
      std::variant<const Framework::Selector::ScalarField<MHDFloat> *, const Framework::Selector::ScalarField<MHDComplex> *>  pField;
      std::variant<std::shared_ptr<Framework::Selector::ScalarField<MHDFloat>>, std::shared_ptr<Framework::Selector::ScalarField<MHDComplex>>>  tField;
      SharedResolution spFieldRes;
      for(auto& fIt: r)
      {
         DebuggerMacro_msg("Add " + ModelOperator::Coordinator::tag(opId) + " term from " + PhysicalNames::Coordinator::tag(fIt.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(fIt.second)) + ")", 7);

         if(fIt.second == FieldComponents::Spectral::SCALAR)
         {
            std::visit(
                  [&](auto&& p)
                  {
                     pField = &p->dom(0).perturbation();
                     spFieldRes = p->dom(0).spRes();
                  }, scalVar.find(fIt.first)->second);
         } else
         {
            std::visit(
                  [&](auto&& p)
                  {
                     pField = &p->dom(0).perturbation().comp(fIt.second);
                     spFieldRes = p->dom(0).spRes();
                  }, vectVar.find(fIt.first)->second);
         }
         if(spEq->res().sim().ss().id() != spFieldRes->sim().ss().id())
         {
            QuICCEnv().synchronize();

            #ifdef QUICC_MPI
               using namespace QuICC::Transpose::Mpi;
            #else
               using namespace QuICC::Transpose::Cpu;
            #endif
            using namespace QuICC::Transpose;

            auto spSetup = spEq->res().spSpectralSetup();
            if(pField.index() == 0)
            {
               auto spF = std::make_shared<Framework::Selector::ScalarField<MHDFloat>>(spSetup);
               tField = spF;
               using VoutTy = View::View<MHDFloat, View::DCCSC3D>;
               using VinTy = View::View<MHDFloat, View::DCCSC3D>;

               // Create transpose operator
               #ifdef QUICC_MPI
                  auto transposeOp = std::make_unique<
                     OpGrouped<std::vector<VoutTy>, std::vector<VinTy>, p021_t>>(spSetup->mem());
               #else
                  auto transposeOp = std::make_unique<
                     OpGrouped<std::vector<VoutTy>, std::vector<VinTy>, p021_t>>();
               #endif

               // Apply transpose
               std::vector<VoutTy> viewsOut = {spF->rGlobalView()};
               std::vector<VinTy> viewsIn = {std::get<0>(pField)->globalView()};
               transposeOp->apply(viewsOut, viewsIn);

               // Update pointer
               pField = spF.get();
            }
            else
            {
               auto spF = std::make_shared<Framework::Selector::ScalarField<MHDComplex>>(spSetup);
               tField = spF;
               using VoutTy = View::View<MHDComplex, View::DCCSC3D>;
               using VinTy = View::View<MHDComplex, View::DCCSC3D>;

               // Create transpose operator
               #ifdef QUICC_MPI
                  auto transposeOp = std::make_unique<
                     OpGrouped<std::vector<VoutTy>, std::vector<VinTy>, p021_t>>(spSetup->mem());
               #else
                  auto transposeOp = std::make_unique<
                     OpGrouped<std::vector<VoutTy>, std::vector<VinTy>, p021_t>>();
               #endif

               // Apply operator
               std::vector<VoutTy> viewsOut = {spF->rGlobalView()};
               std::vector<VinTy> viewsIn = {std::get<1>(pField)->globalView()};
               transposeOp->apply(viewsOut, viewsIn);

               // Update pointer
               pField = spF.get();
            }

            QuICCEnv().synchronize();
         }

         // Get explicit input
         std::visit(
               [&](auto&& p)
               {
                  for(std::size_t i = 0; i < (*solveIt)->nSystem(); i++)
                  {
                     addExp(spEq, fIt, i, *p);
                  }
               }, pField);
      }
   }

   template <typename TSolverIt,typename TEq> void computeSolverInput(const TEq spEq, const SpectralFieldId id, const TSolverIt solveIt, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>&, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>&)
   {
      const auto& cinfo = spEq->couplingInfo(id.second);
      const bool hasQI  = (cinfo.hasNonlinear() && cinfo.hasQuasiInverse());
      const auto compId = id.second;

      // Get timestep input
      for(std::size_t i = 0; i < (*solveIt)->nSystem(); i++)
      {
         if(hasQI)
         {
            if(spEq->hasQID(compId))
            {
               const SparseMatrix& op = spEq->template quasiInverse<SparseMatrix>(compId, i);
               std::visit(
                     [&](auto&& p)
                     {
                        Equations::copyNonlinear(spEq->res(), cinfo, p->dom(0).perturbation().comp(id.second), op, (*solveIt)->rRHSData(i), i, (*solveIt)->startRow(id,i));
                     }, spEq->spUnknown());
            }
            else if(spEq->hasQIZ(compId))
            {
               const SparseMatrixZ& op = spEq->template quasiInverse<SparseMatrixZ>(compId, i);
               std::visit(
                     [&](auto&& p)
                     {
                        Equations::copyNonlinear(spEq->res(), cinfo, p->dom(0).perturbation().comp(id.second), op, (*solveIt)->rRHSData(i), i, (*solveIt)->startRow(id,i));
                     }, spEq->spUnknown());
            }
         }
         else
         {
            // Copy field values into solver input
            std::visit(
                  [&](auto&& p)
                  {
                     Equations::copyNonlinear(spEq->res(), cinfo, p->dom(0).perturbation().comp(id.second), (*solveIt)->rRHSData(i), i, (*solveIt)->startRow(id,i));
                  }, spEq->spUnknown());
         }

         // Add source term
         if(cinfo.hasSource())
         {
            std::visit(
                  [&](auto&& p)
                  {
                     Equations::addSource(spEq->res(), cinfo, spEq->spSourceKernel(id.second), p->dom(0).perturbation().comp(id.second), (*solveIt)->rRHSData(i), i, (*solveIt)->startRow(id,i));
                  }, spEq->spUnknown());
         }

         // If required set inhomogenous boundary condition value
         if(cinfo.hasBoundaryValue())
         {
            // Set boundary value
            std::visit(
                  [&](auto&& p)
                  {
                     Equations::setBoundaryValue(spEq->res(), cinfo, spEq->spBoundaryKernel(id.second), p->dom(0).perturbation().comp(id.second), (*solveIt)->rInhomogeneous(i), i, (*solveIt)->startRow(id,i));
                  }, spEq->spUnknown());
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

} // Solver
} // QuICC

#endif // QUICC_SOLVER_SPARSECOORDINATORDATA_HPP
