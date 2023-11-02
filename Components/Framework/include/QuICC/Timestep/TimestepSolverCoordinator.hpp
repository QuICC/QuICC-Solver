/**
 * @file TimestepSolverCoordinator.hpp
 * @brief Implementation of a interface to the RK CB schemes
 */

#ifndef QUICC_TIMESTEP_TIMESTEPSOLVERCOORDINATOR_HPP
#define QUICC_TIMESTEP_TIMESTEPSOLVERCOORDINATOR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "Types/Math.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
#include "QuICC/ModelOperator/SplitBoundaryValue.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/SparseSolvers/SparseLinearCoordinatorBase.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Implementation of general timestepper structure
    */
   template <template <class,class,template <class> class> class TSolver> class TimestepSolverCoordinator: public Solver::SparseLinearCoordinatorBase<TSolver>
   {
      public:
         /**
          * @brief Constructor
          */
         TimestepSolverCoordinator() = default;

         /**
          * @brief Destructor
          */
         virtual ~TimestepSolverCoordinator() = default;

         /**
          * @brief Initialise solver coordinator
          *
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         template <typename TScheme> void init(const MHDFloat dt, const typename TimestepSolverCoordinator<TSolver>::ScalarEquation_range& scalEq, const typename TimestepSolverCoordinator<TSolver>::VectorEquation_range& vectEq, std::shared_ptr<TScheme> spScheme);

         /**
          * @brief Update timestep
          */
         void updateTimestep(const MHDFloat dt);

      protected:
         /**
          * @brief Build the real operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(typename TimestepSolverCoordinator<TSolver>::SharedRealSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Build the complex operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(typename TimestepSolverCoordinator<TSolver>::SharedComplexSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

      private:
         using Solver::SparseLinearCoordinatorBase<TSolver>::init;

         /**
          * @brief Current timestep
          */
         MHDFloat mDt;
   };

   /**
    * @brief Small wrapper for a generic implementation of the solver matrix construction
    */
   template <typename TStepper> void buildTimestepMatrixWrapper(typename std::shared_ptr<TStepper > spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const MHDFloat dt, const int idx);

   template <typename TStepper> void buildTimestepMatrixWrapper(typename std::shared_ptr<TStepper > spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const MHDFloat dt, const int idx)
   {
      std::map<std::size_t, DecoupledZSparse> ops;

      bool isSplit = spEq->couplingInfo(comp).isSplitEquation();

      // Compute model's linear operator (without Tau lines)
      ops.insert(std::make_pair(ModelOperator::ImplicitLinear::id(), DecoupledZSparse()));
      spEq->buildModelMatrix(ops.find(ModelOperator::ImplicitLinear::id())->second, ModelOperator::ImplicitLinear::id(), comp, idx, ModelOperatorBoundary::SolverNoTau::id());
      // Compute model's time operator (without Tau lines)
      ops.insert(std::make_pair(ModelOperator::Time::id(), DecoupledZSparse()));
      spEq->buildModelMatrix(ops.find(ModelOperator::Time::id())->second, ModelOperator::Time::id(), comp, idx, ModelOperatorBoundary::SolverNoTau::id());
      // Compute model's tau line boundary operator
      ops.insert(std::make_pair(ModelOperator::Boundary::id(), DecoupledZSparse()));
      spEq->buildModelMatrix(ops.find(ModelOperator::Boundary::id())->second, ModelOperator::Boundary::id(), comp, idx, ModelOperatorBoundary::SolverHasBc::id());

      // If equation was split into two lower order systems
      if(isSplit)
      {
         // Compute model's split linear operator (without Tau lines)
         auto id = ModelOperator::SplitImplicitLinear::id();
         ops.insert(std::make_pair(id, DecoupledZSparse()));
         spEq->buildModelMatrix(ops.find(id)->second, id, comp, idx, ModelOperatorBoundary::SolverNoTau::id());

         // Compute model's tau line boundary operator for split operator
         id = ModelOperator::SplitBoundary::id();
         ops.insert(std::make_pair(id, DecoupledZSparse()));
         spEq->buildModelMatrix(ops.find(id)->second, id, comp, idx, ModelOperatorBoundary::SolverHasBc::id());

         // Compute model's tau line boundary value for split operator
         id = ModelOperator::SplitBoundaryValue::id();
         ops.insert(std::make_pair(id, DecoupledZSparse()));
         spEq->buildModelMatrix(ops.find(id)->second, id, comp, idx, ModelOperatorBoundary::SolverHasBc::id());
      }

      // Let the timestepper build the right operators
      spSolver->buildOperators(idx, ops, dt, spEq->couplingInfo(comp).systemN(idx));

      // Solver is initialized
      spSolver->setInitialized();
   }

   template <template <class,class,template <class> class> class TSolver> void TimestepSolverCoordinator<TSolver>::buildSolverMatrix(typename TimestepSolverCoordinator<TSolver>::SharedRealSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      buildTimestepMatrixWrapper(spSolver, spEq, comp, this->mDt, idx);
   }

   template <template <class,class,template <class> class> class TSolver> void TimestepSolverCoordinator<TSolver>::buildSolverMatrix(typename TimestepSolverCoordinator<TSolver>::SharedComplexSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      buildTimestepMatrixWrapper(spSolver, spEq, comp, this->mDt, idx);
   }

   template <template <class,class,template <class> class> class TSolver> void TimestepSolverCoordinator<TSolver>::updateTimestep(const MHDFloat dt)
   {
      this->mDt = dt;
   }

   template <template <class,class,template <class> class> class TSolver> template <typename TScheme> void TimestepSolverCoordinator<TSolver>::init(const MHDFloat dt, const typename TimestepSolverCoordinator<TSolver>::ScalarEquation_range& scalEq, const typename TimestepSolverCoordinator<TSolver>::VectorEquation_range& vectEq, std::shared_ptr<TScheme> spScheme)
   {
      this->mDt = dt;

      // Create solvers
      this->createCoordSolvers(scalEq, vectEq);

      // Set Timestepper scheme
      for(auto solIt = this->mRealSolvers.begin(); solIt != this->mRealSolvers.end(); ++ solIt)
      {
         (*solIt)->setScheme(spScheme);
      }

      for(auto solIt = this->mComplexSolvers.begin(); solIt != this->mComplexSolvers.end(); ++ solIt)
      {
         (*solIt)->setScheme(spScheme);
      }

      // Initialise solver storage
      this->initCoordSolvers(scalEq, vectEq);
   }

} // Timestep
} // QuICC

#endif // QUICC_TIMESTEP_TIMESTEPSOLVERCOORDINATOR_HPP
