/**
 * @file SparseCoordinatorBase.hpp
 * @brief Implementation of the base for a general sparse solver coordinator
 */

#ifndef QUICC_SOLVER_SPARSECOORDINATORBASE_HPP
#define QUICC_SOLVER_SPARSECOORDINATORBASE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/IteratorRange.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/SparseSolvers/SparseCoordinatorData.hpp"

namespace QuICC {

namespace Solver {

   /**
    * @brief Implementation of the base for a general sparse solver coordinator
    */
   template <template <class,class,template <typename> class> class TSolver> class SparseCoordinatorBase: public SparseCoordinatorData<TSolver>
   {
      public:
         /// Typedef for a shared scalar equation iterator
         typedef typename std::vector<Equations::SharedIScalarEquation>::iterator   ScalarEquation_iterator;

         /// Typedef for a shared vector equation iterator
         typedef typename std::vector<Equations::SharedIVectorEquation>::iterator   VectorEquation_iterator;

         /// Typedef for a shared scalar equation range
         typedef std::pair<ScalarEquation_iterator, ScalarEquation_iterator>  ScalarEquation_range;

         /// Typedef for a shared vector equation range
         typedef std::pair<VectorEquation_iterator, VectorEquation_iterator>  VectorEquation_range;

         /**
          * @brief Constructor
          */
         SparseCoordinatorBase();

         /**
          * @brief Destructor
          */
         virtual ~SparseCoordinatorBase();

         /**
          * @brief Initialise solver coordinator
          *
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         virtual void init(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq) = 0;

         /**
          * @brief Update equation explicit linear input to solver
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void getExplicitInput(const std::size_t opId, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const typename SparseCoordinatorBase<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorBase<TSolver>::VectorVariable_map& vectVar);

         /**
          * @brief Update equation input to solver
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const typename SparseCoordinatorBase<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorBase<TSolver>::VectorVariable_map& vectVar);

         /**
          * @brief Update equation unkowns with solver output
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          */
         void transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Get error measure
          */
         MHDFloat error() const;

         /**
          * @brief Check if step is done
          */
         bool finishedStep() const;

      protected:
         /**
          * @brief Create the correct solver
          */
         void createSolver(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Create storage for the solvers
          */
         void createStorage(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Initialise the solution
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void initSolution(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Flag to signal end of computation
          */
         bool mFinished;

         /**
          * @brief Computation error
          */
         MHDFloat mError;

      private:
   };

   template <template <class,class,template <class> class> class TSolver> SparseCoordinatorBase<TSolver>::SparseCoordinatorBase()
      : SparseCoordinatorData<TSolver>(), mFinished(false), mError(-1.0)
   {
   }

   template <template <class,class,template <class> class> class TSolver> SparseCoordinatorBase<TSolver>::~SparseCoordinatorBase()
   {
   }

   template <template <class,class,template <class> class> class TSolver> MHDFloat SparseCoordinatorBase<TSolver>::error() const
   {
      return this->mError;
   }

   template <template <class,class,template <class> class> class TSolver> bool SparseCoordinatorBase<TSolver>::finishedStep() const
   {
      return this->mFinished;
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorBase<TSolver>::createSolver(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp)
   {
      DebuggerMacro_msg("Creating solver for " + PhysicalNames::Coordinator::tag(spEq->name()) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(comp)) + ")", 2);

      // System has a complex operator
      if(spEq->couplingInfo(comp).isComplex())
      {
         typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator solIt;
         this->addSolver(solIt, spEq->couplingInfo(comp).solverIndex(), spEq->couplingInfo(comp).fieldStart(), spEq->solveTiming());

      // System has a real operator
      } else
      {
         typename SparseCoordinatorBase<TSolver>::RealSolver_iterator solIt;
         this->addSolver(solIt, spEq->couplingInfo(comp).solverIndex(), spEq->couplingInfo(comp).fieldStart(), spEq->solveTiming());
      }

      DebuggerMacro_msg("... done", 2);
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorBase<TSolver>::createStorage(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp)
   {
      DebuggerMacro_msg("Creating storage for solver for " + PhysicalNames::Coordinator::tag(spEq->name()) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(comp)) + ")", 2);

      // ID of the current field
      SpectralFieldId myId = std::make_pair(spEq->name(),comp);

      // Index of solver
      int myIdx = spEq->couplingInfo(myId.second).solverIndex();

      // System has a complex operator
      if(spEq->couplingInfo(myId.second).isComplex())
      {
         setupSolverStorage<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, spEq, myIdx, myId);

      // System has a real operator
      } else
      {
         setupSolverStorage<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, spEq, myIdx, myId);
      }

      DebuggerMacro_msg("... done", 2);
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorBase<TSolver>::transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Storage for identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      for(auto& scalEqIt: make_range(scalEq))
      {
         if(scalEqIt->solveTiming() == this->solveTime())
         {
            // Get field identity
            myId = std::make_pair(scalEqIt->name(), FieldComponents::Spectral::SCALAR);

            // Get index of solver
            int myIdx = scalEqIt->couplingInfo(myId.second).solverIndex();

            // System operator is complex
            if(scalEqIt->couplingInfo(myId.second).isComplex())
            {
               storeSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, scalEqIt, myIdx, myId);

               // System operator is real real
            } else
            {
               storeSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, scalEqIt, myIdx, myId);
            }
         }
      }

      // Loop over all vector equations
      for(auto& vectEqIt: make_range(vectEq))
      {
         if(vectEqIt->solveTiming() == this->solveTime())
         {
            // Loop over the vector equation components
            for(auto& compIt: make_range(vectEqIt->spectralRange()))
            {
               // Get field identity
               myId = std::make_pair(vectEqIt->name(), compIt);

               // Get index of solver
               int myIdx = vectEqIt->couplingInfo(myId.second).solverIndex();

               // System operator is complex
               if(vectEqIt->couplingInfo(myId.second).isComplex())
               {
                  storeSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, vectEqIt, myIdx, myId);

               // System operator is real
               } else
               {
                  storeSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, vectEqIt, myIdx, myId);
               }
            }
         }
      }
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorBase<TSolver>::initSolution(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      for(auto& scalEqIt: make_range(scalEq))
      {
         // Get field identity
         myId = std::make_pair(scalEqIt->name(), FieldComponents::Spectral::SCALAR);

         // Get index of solver
         int myIdx = scalEqIt->couplingInfo(myId.second).solverIndex();

         // System operator is complex
         if(scalEqIt->couplingInfo(myId.second).isComplex())
         {
            initSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, scalEqIt, myIdx, myId);

         // System operator is real
         } else
         {
            initSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, scalEqIt, myIdx, myId);
         }
      }

      // Loop over all vector equations
      for(auto& vectEqIt: make_range(vectEq))
      {
         for(auto& compIt: make_range(vectEqIt->spectralRange()))
         {
            // Get field identity
            myId = std::make_pair(vectEqIt->name(), compIt);

            // Get index of solver
            int myIdx = vectEqIt->couplingInfo(myId.second).solverIndex();

            // System operator is complex
            if(vectEqIt->couplingInfo(myId.second).isComplex())
            {
               initSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, vectEqIt, myIdx, myId);

            // System operator is real
            } else
            {
               initSolverSolution<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, vectEqIt, myIdx, myId);
            }
         }
      }
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorBase<TSolver>::getExplicitInput(const std::size_t opId, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const typename SparseCoordinatorBase<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorBase<TSolver>::VectorVariable_map& vectVar)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      for(auto& scalEqIt: make_range(scalEq))
      {
         // Get field identity
         myId = std::make_pair(scalEqIt->name(), FieldComponents::Spectral::SCALAR);

         // Get index of solver
         int myIdx = scalEqIt->couplingInfo(myId.second).solverIndex();

         // System operator is complex
         if(scalEqIt->couplingInfo(myId.second).isComplex())
         {
            getExplicitSolverInput<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, scalEqIt, myIdx, myId, opId, scalVar, vectVar);

            // System operator is real
         } else
         {
            getExplicitSolverInput<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, scalEqIt, myIdx, myId, opId, scalVar, vectVar);
         }
      }

      // Loop over all vector equations
      for(auto& vectEqIt: make_range(vectEq))
      {
         for(auto& compIt: make_range(vectEqIt->spectralRange()))
         {
            // Get field identity for first component
            myId = std::make_pair(vectEqIt->name(), compIt);

            // Get index of solver
            int myIdx = vectEqIt->couplingInfo(myId.second).solverIndex();

            // Linear solve matrices are complex
            if(vectEqIt->couplingInfo(myId.second).isComplex())
            {
               getExplicitSolverInput<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, vectEqIt, myIdx, myId, opId, scalVar, vectVar);

               // Linear solve matrices are real
            } else
            {
               getExplicitSolverInput<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, vectEqIt, myIdx, myId, opId, scalVar, vectVar);
            }
         }
      }
   }

   template <template <class,class,template <class> class> class TSolver> void SparseCoordinatorBase<TSolver>::getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const typename SparseCoordinatorBase<TSolver>::ScalarVariable_map& scalVar, const typename SparseCoordinatorBase<TSolver>::VectorVariable_map& vectVar)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      for(auto& scalEqIt: make_range(scalEq))
      {
         if(scalEqIt->solveTiming() == this->solveTime())
         {
            // Get field identity
            myId = std::make_pair(scalEqIt->name(), FieldComponents::Spectral::SCALAR);

            // Get index of solver
            int myIdx = scalEqIt->couplingInfo(myId.second).solverIndex();

            // System operator is complex
            if(scalEqIt->couplingInfo(myId.second).isComplex())
            {
               getSolverInput<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, scalEqIt, myIdx, myId, scalVar, vectVar);

            // System operator is real
            } else
            {
               getSolverInput<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, scalEqIt, myIdx, myId, scalVar, vectVar);
            }
         }
      }

      // Loop over all vector equations
      for(auto& vectEqIt: make_range(vectEq))
      {
         if(vectEqIt->solveTiming() == this->solveTime())
         {
            for(auto& compIt: make_range(vectEqIt->spectralRange()))
            {
               // Get field identity for first component
               myId = std::make_pair(vectEqIt->name(), compIt);

               // Get index of solver
               int myIdx = vectEqIt->couplingInfo(myId.second).solverIndex();

               // Linear solve matrices are complex
               if(vectEqIt->couplingInfo(myId.second).isComplex())
               {
                  getSolverInput<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, vectEqIt, myIdx, myId, scalVar, vectVar);

                  // Linear solve matrices are real
               } else
               {
                  getSolverInput<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, vectEqIt, myIdx, myId, scalVar, vectVar);
               }
            }
         }
      }
   }
} // Solver
} // QuICC

#endif // QUICC_SOLVER_SPARSECOORDINATORBASE_HPP
