/**
 * @file SparseLinearCoordinatorBase.hpp
 * @brief Implementation of the base for a general sparse linear solver coordinator
 */

#ifndef QUICC_SOLVER_SPARSELINEARCOORDINATORBASE_HPP
#define QUICC_SOLVER_SPARSELINEARCOORDINATORBASE_HPP

// System includes
//
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// Project includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/SparseSolvers/SparseCoordinatorBase.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/Timers/StageTimer.hpp"

namespace QuICC {

namespace Solver {

   /**
    * @brief Implementation of the base for a general sparse linear solver coordinator
    */
   template <template <class,class,template <class> class> class TSolver> class SparseLinearCoordinatorBase: public SparseCoordinatorBase<TSolver>
   {
      public:
         /// Typedef for a shared scalar equation range
         typedef typename SparseCoordinatorBase<TSolver>::ScalarEquation_range  ScalarEquation_range;

         /// Typedef for a shared vector equation range
         typedef typename SparseCoordinatorBase<TSolver>::VectorEquation_range  VectorEquation_range;

         /// Typedef for a shared real operator solver
         typedef typename SparseCoordinatorBase<TSolver>::SharedRealSolverType  SharedRealSolverType;

         /// Typedef for a shared complex operator solver
         typedef typename SparseCoordinatorBase<TSolver>::SharedComplexSolverType  SharedComplexSolverType;

         /**
          * @brief Constructor
          */
         SparseLinearCoordinatorBase();

         /**
          * @brief Destructor
          */
         virtual ~SparseLinearCoordinatorBase() = default;

         /**
          * @brief Initialise solver coordinator
          *
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void init(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Build the solver matrices independently of solver type
          *
          * \mhdBug Should not be public
          */
         template <typename TSolverIt> void buildSolverMatrices(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt);

         /**
          * @brief Solve all the linear systems
          */
         void solveSystems();

      protected:
         /**
          * @brief create solvers
          *
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void createCoordSolvers(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Initialise solvers
          *
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void initCoordSolvers(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Compute (coupled) matrices
          */
         void createMatrices(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Build the real operator
          *
          * @param spSolver   Shared sparse real solver
          * @param spEq    Shared pointer to equation
          * @param comp    Field component
          * @param idx     Matrix index
          */
         virtual void buildSolverMatrix(SharedRealSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx) = 0;

         /**
          * @brief Build the complex operator
          *
          * @param spSolver   Shared sparse real solver
          * @param spEq    Shared pointer to equation
          * @param comp    Field component
          * @param idx     Matrix index
          */
         virtual void buildSolverMatrix(SharedComplexSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx) = 0;

         /**
          * @brief Factorization time
          */
         MHDFloat mFactorTime;

      private:
   };

   template <template <class,class,template <class> class> class TSolver, typename TSolverIt> void buildSolverMatricesSolver(SparseLinearCoordinatorBase<TSolver>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId id);

   //
   //
   //

   template <template <class,class,template <class> class> class TSolver> template <typename TSolverIt> void SparseLinearCoordinatorBase<TSolver>::buildSolverMatrices(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt)
   {
      // Number of linear systems
      int nSystems = spEq->couplingInfo(id.second).nSystems();

      // Reserve storage for matrices and initialise vectors
      (*solveIt)->initMatrices(nSystems);

      // Build the solver matrices
      for(int i = 0; i < nSystems; i++)
      {
         // Build LHS solver matrix
         this->buildSolverMatrix(*solveIt, spEq, id.second, i);
      }
   }

   template <template <class,class,template <class> class> class TSolver> SparseLinearCoordinatorBase<TSolver>::SparseLinearCoordinatorBase()
      : SparseCoordinatorBase<TSolver>(), mFactorTime(-1.0)
   {
   }

   template <template <class,class,template <class> class> class TSolver> void SparseLinearCoordinatorBase<TSolver>::createCoordSolvers(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      StageTimer  stage;
      //
      // Create real/complex solvers
      //

      stage.start("creating linear solvers", 1);

      // Loop over all scalar equations
      for(auto& scalEqIt: make_range(scalEq))
      {
         // Get type information for the linear solvers
         this->createSolver(scalEqIt, FieldComponents::Spectral::SCALAR);
      }

      // Loop over all vector equations
      for(auto& vectEqIt: make_range(vectEq))
      {
         // Get type information for the linear solvers
         for(auto& compIt: make_range(vectEqIt->spectralRange()))
         {
            this->createSolver(vectEqIt, compIt);
         }
      }

      stage.done();
   }

   template <template <class,class,template <class> class> class TSolver> void SparseLinearCoordinatorBase<TSolver>::initCoordSolvers(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      StageTimer  stage;

      //
      // Initialise the solver storage
      //

      stage.start("initializing solver storage", 1);
      // Loop over all scalar equations
      for(auto& scalEqIt: make_range(scalEq))
      {
         // Create storage
         this->createStorage(scalEqIt, FieldComponents::Spectral::SCALAR);
      }

      // Loop over all vector equations
      for(auto& vectEqIt: make_range(vectEq))
      {
         // Create storage
         for(auto& compIt: make_range(vectEqIt->spectralRange()))
         {
            this->createStorage(vectEqIt, compIt);
         }
      }

      // Initialise the start rows
      this->initStartRow();

      stage.done();

      //
      // Create the problem matrices
      //

      stage.start("building solver matrices", 1);

      // Loop over all scalar equations
      for(auto& scalEqIt: make_range(scalEq))
      {
         // Create (coupled) matrices
         this->createMatrices(scalEqIt, FieldComponents::Spectral::SCALAR);
      }

      // Loop over all vector equations
      for(auto& vectEqIt: make_range(vectEq))
      {
         // Create (coupled) matrices
         for(auto& compIt: make_range(vectEqIt->spectralRange()))
         {
            this->createMatrices(vectEqIt, compIt);
         }
      }

      stage.done();

      //
      // Initialise the solvers and the initial state
      //

      stage.start("factorizing solver matrices", 1);
      TimerMacro timer(true);

      // Initialise solvers from complex equation steppers
      initSolvers<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this);

      // Initialise solvers from real equation steppers
      initSolvers<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this);

      timer.stop();
      this->mFactorTime = timer.time();
      stage.done();

      stage.start("initializing solutions", 1);

      // Initialise with initial state
      this->initSolution(scalEq, vectEq);

      stage.done();
   }

   template <template <class,class,template <class> class> class TSolver> void SparseLinearCoordinatorBase<TSolver>::init(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Create the solvers
      this->createCoordSolvers(scalEq, vectEq);

      // Initialize the solvers
      this->initCoordSolvers(scalEq, vectEq);
   }

   template <template <class,class,template <class> class> class TSolver> void SparseLinearCoordinatorBase<TSolver>::createMatrices(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp)
   {
      DebuggerMacro_msg("creating operators for " + PhysicalNames::Coordinator::tag(spEq->name()) + "(" + std::to_string(comp) + ")", 2);

      // ID of the current field
      SpectralFieldId myId = std::make_pair(spEq->name(),comp);

      // Get solver index
      int myIdx = spEq->couplingInfo(myId.second).solverIndex();

      // System operator is complex
      if(spEq->couplingInfo(myId.second).isComplex())
      {
         buildSolverMatricesSolver<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, spEq, myIdx, myId);

      // System operator is real
      } else
      {
         buildSolverMatricesSolver<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, spEq, myIdx, myId);
      }

      DebuggerMacro_msg("... done", 2);
   }

   template <template <class,class,template <class> class> class TSolver> void SparseLinearCoordinatorBase<TSolver>::solveSystems()
   {
      // Solve complex operator, complex field linear systems
      std::pair<bool,MHDFloat> zStatus = solveSolvers<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this);

      // Solve real operator, complex field linear systems
      std::pair<bool,MHDFloat> dStatus = solveSolvers<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this);

      this->mFinished = zStatus.first || dStatus.first;

      if(this->mFinished)
      {
         this->mError = std::max(zStatus.second, dStatus.second);
      }
   }

   //
   //
   //

   template <template <class,class,template <class> class> class TSolver, typename TSolverIt> void buildSolverMatricesSolver(SparseLinearCoordinatorBase<TSolver>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId id)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt, idx);

      // Build solver matrices
      if(! (*solIt)->isInitialized())
      {
         coord.buildSolverMatrices(spEq, id, solIt);
      }
   }
} // Solver
} // QuICC

#endif // QUICC_SOLVER_SPARSELINEARCOORDINATORBASE_HPP
