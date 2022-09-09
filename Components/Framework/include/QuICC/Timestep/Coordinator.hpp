/**
 * @file Coordinator.hpp
 * @brief Implementation of a general timestep coordinator structure
 */

#ifndef QUICC_TIMESTEP_COORDINATOR_HPP
#define QUICC_TIMESTEP_COORDINATOR_HPP

// Configuration includes
//

// Debug includes
//
#include "QuICC/Debug/DebuggerMacro.h"

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/SparseSolvers/SparseLinearCoordinatorBase.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/Timestep/IScheme.hpp"
#include "QuICC/TypeSelectors/TimeSchemeSelector.hpp"
#include "QuICC/Io/Ascii/CflWriter.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Implementation of general timestepper structure
    */
   class Coordinator: public Solver::SparseLinearCoordinatorBase<TimeSchemeSelector::ImplementationType>
   {
      public:
         /// Typedef for parent coordinator
         typedef typename Solver::SparseLinearCoordinatorBase<TimeSchemeSelector::ImplementationType> ParentCoordinator;

         /// Typedef for a shared scalar equation range
         typedef typename ParentCoordinator::ScalarEquation_range  ScalarEquation_range;

         /// Typedef for a shared vector equation range
         typedef typename ParentCoordinator::VectorEquation_range  VectorEquation_range;

         /// Typedef for a shared scalar variable map
         typedef typename ParentCoordinator::ScalarVariable_map  ScalarVariable_map;

         /// Typedef for a shared vector variable map
         typedef typename ParentCoordinator::VectorVariable_map  VectorVariable_map;

         /// Typedef for a shared real operator solver
         typedef typename ParentCoordinator::SharedRealSolverType  SharedRealSolverType;

         /// Typedef for a shared complex operator solver
         typedef typename ParentCoordinator::SharedComplexSolverType  SharedComplexSolverType;

         /**
          * @brief Constructor
          */
         Coordinator();

         /**
          * @brief Destructor
          */
         ~Coordinator();

         /**
          * @brief Initialise timestepper
          *
          * @param time    Initial time value
          * @param cfl     Initial CFL timestep
          * @param error   Max error allowed during timestep
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void init(const MHDFloat time, const Matrix& cfl, const MHDFloat maxError, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Tune adaptive timestepper
          *
          * \mhdBug Not fully implemented
          */
         void tuneAdaptive(const MHDFloat time);

         /**
          * @brief Adapt the timestep used
          *
          * @param cfl     CFL conditions
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void adaptTimestep(const Matrix& cfl, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Update control status
          */
         void update();

         /**
          * @brief Compute (partial) forward step
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          * @param scalVar Shared scalar variables
          * @param vectVar Shared vector variables
          */
         void stepForward(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar);

         /**
          * @brief Get current simulation time
          */
         MHDFloat time() const;

         /**
          * @brief Get current simulation timestep
          */
         MHDFloat timestep() const;

         /**
          * @brief Print timestepper information to stream
          *
          * @param stream  Output stream
          */
         void printInfo(std::ostream& stream);

      protected:
         /**
          * @brief Build the real operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(Coordinator::SharedRealSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Build the complex operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(Coordinator::SharedComplexSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

      private:
         using ParentCoordinator::init;

         /**
          * @brief Update time dependence
          */
         void updateMatrices();

         /**
          * @brief Minimum number of constant timestep before step size increase
          */
         const MHDFloat mcMinCnst;

         /**
          * @brief Maximum timestep jump per step (See Soederlind)
          */
         const MHDFloat mcMaxJump;

         /**
          * @brief No update window for timestep increase
          */
         const MHDFloat mcUpWindow;

         /**
          * @brief Minimal timestep allowed before simulation abort
          */
         const MHDFloat mcMinDt;

         /**
          * @brief Maximum timestep allowed
          */
         const MHDFloat mcMaxDt;

         /**
          * @brief Maximum error allowed
          */
         MHDFloat mMaxError;

         /**
          * @brief Previous timestep length
          */
         MHDFloat mOldDt;

         /**
          * @brief Timestep length
          */
         Matrix mDt;

         /**
          * @brief Current time
          */
         MHDFloat mTime;

         /**
          * @brief Current reference time
          */
         MHDFloat mRefTime;

         /**
          * @brief Constant timestep steps
          */
         MHDFloat mCnstSteps;

         /**
          * @brief Constant timestep steps
          */
         MHDFloat mStepTime;

         /**
          * @brief Shared CFL writer
          */
         Io::Ascii::SharedCflWriter   mspIo;

         /**
          * @brief Interface to timestepping scheme
          */
         SharedIScheme  mspScheme;
   };

   /**
    * @brief Small wrapper for a generic implementation of the solver matrix construction
    */
   template <typename TStepper> void buildTimestepMatrixWrapper(typename std::shared_ptr<TStepper > spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const MHDFloat dt, const int idx);

   template <typename TStepper> void buildTimestepMatrixWrapper(typename std::shared_ptr<TStepper > spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const MHDFloat dt, const int idx)
   {
      std::map<std::size_t, DecoupledZSparse> ops;

      // Compute model's linear operator (without Tau lines)
      ops.insert(std::make_pair(ModelOperator::ImplicitLinear::id(), DecoupledZSparse()));
      spEq->buildModelMatrix(ops.find(ModelOperator::ImplicitLinear::id())->second, ModelOperator::ImplicitLinear::id(), comp, idx, ModelOperatorBoundary::SolverNoTau::id());
      // Compute model's time operator (without Tau lines)
      ops.insert(std::make_pair(ModelOperator::Time::id(), DecoupledZSparse()));
      spEq->buildModelMatrix(ops.find(ModelOperator::Time::id())->second, ModelOperator::Time::id(), comp, idx, ModelOperatorBoundary::SolverNoTau::id());
      // Compute model's tau line boundary operator
      ops.insert(std::make_pair(ModelOperator::Boundary::id(), DecoupledZSparse()));
      spEq->buildModelMatrix(ops.find(ModelOperator::Boundary::id())->second, ModelOperator::Boundary::id(), comp, idx, ModelOperatorBoundary::SolverHasBc::id());

      // Let the timestepper build the right operators
      spSolver->buildOperators(idx, ops, dt, spEq->couplingInfo(comp).systemN(idx));

      // Solver is initialized
      spSolver->setInitialized();
   }

}
}

#endif // QUICC_TIMESTEP_COORDINATOR_HPP
