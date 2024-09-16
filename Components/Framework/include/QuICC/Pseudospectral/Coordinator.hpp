/**
 * @file Coordinator.hpp
 * @brief High level pseudospectral coordinator
 */

#ifndef QUICC_PSEUDOSPECTRAL_COORDINATOR_HPP
#define QUICC_PSEUDOSPECTRAL_COORDINATOR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/LoadSplitter/Algorithms/SplittingDescription.hpp"
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "QuICC/SpatialScheme/IBuilder.hpp"
#include "QuICC/Equations/EquationParameters.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/SparseSolvers/SparseTrivialCoordinator.hpp"
#include "QuICC/SparseSolvers/SparseLinearCoordinator.hpp"
#include "QuICC/Io/Config/ConfigurationReader.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"
#include "QuICC/TypeSelectors/ParallelSelector.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"
#include "QuICC/TransformGroupers/IForwardGrouper.hpp"
#include "QuICC/TransformGroupers/IBackwardGrouper.hpp"
#include "QuICC/Timestep/Coordinator.hpp"
#include "QuICC/Diagnostics/Coordinator.hpp"
#include "Graph/Jit.hpp"

namespace QuICC {

namespace Pseudospectral {

   /**
    * @brief High level pseudospectral coordinator
    */
   class Coordinator
   {
      public:
         /// Typedef for a IScalarEquation
         typedef Equations::IScalarEquation IScalarEquation;

         /// Typedef for a IVectorEquation
         typedef Equations::IVectorEquation IVectorEquation;

         /// Typedef for a IScalarEquation
         typedef std::shared_ptr<IScalarEquation> SharedIScalarEquation;

         /// Typedef for a IVectorEquation
         typedef std::shared_ptr<IVectorEquation> SharedIVectorEquation;

         /// Typedef for a shared scalar equation iterator
         typedef std::vector<SharedIScalarEquation>::iterator   ScalarEquation_iterator;

         /// Typedef for a shared vector equation iterator
         typedef std::vector<SharedIVectorEquation>::iterator   VectorEquation_iterator;

         /// Typedef for a shared scalar equation range
         typedef std::pair<ScalarEquation_iterator, ScalarEquation_iterator>  ScalarEquation_range;

         /// Typedef for a shared vector equation range
         typedef std::pair<VectorEquation_iterator, VectorEquation_iterator>  VectorEquation_range;

         /// Typedef for a scalar variable
         typedef Framework::Selector::VariantSharedScalarVariable  ScalarVariable;

         /// Typedef for a vector variable
         typedef Framework::Selector::VariantSharedVectorVariable  VectorVariable;

         /// Typedef for a shared scalar variable map
         typedef std::map<std::size_t, ScalarVariable>  ScalarVariable_map;

         /// Typedef for a shared vector variable map
         typedef std::map<std::size_t, VectorVariable>  VectorVariable_map;

         /// Typedef for a shared scalar variable iterator
         typedef ScalarVariable_map::iterator  ScalarVariable_iterator;

         /// Typedef for a shared vector variable iterator
         typedef VectorVariable_map::iterator  VectorVariable_iterator;

         /// Typedef for a shared scalar variable range
         typedef std::pair<ScalarVariable_iterator, ScalarVariable_iterator>  ScalarVariable_range;

         /// Typedef for a shared vector variable range
         typedef std::pair<VectorVariable_iterator, VectorVariable_iterator>  VectorVariable_range;

         /// Typedef for trivial solver coordinator
         typedef Solver::SparseTrivialCoordinator SparseTrivialCoordinator;

         /// Typedef for linear solver coordinator
         typedef Solver::SparseLinearCoordinator SparseLinearCoordinator;

         /**
          * @brief Constructor
          */
         Coordinator() = default;

         /**
          * @brief Destructor
          */
         virtual ~Coordinator() = default;

         /**
          * @brief Get iteration indexes
          */
         const std::set<int>& it() const;

         /**
          * @brief Get integration time
          */
         MHDFloat time() const;

         /**
          * @brief Get integration timestep
          */
         MHDFloat timestep() const;

         /**
          * @brief Get start time
          */
         MHDFloat startTime() const;

         /**
          * @brief Get start timestep
          */
         MHDFloat startTimestep() const;

         /**
          * @brief Print pseudospectral information to stream
          *
          * @param stream  Output stream
          */
         virtual void printInfo(std::ostream& stream);

         /**
          * @brief Evolve pseudospectral equations
          */
         virtual void evolve();

         /**
          * @brief Evolve pseudospectral equations
          *
          * @param it Current iteration
          */
         virtual void evolveUntilPrognostic(const bool finishedStep);

         /**
          * @brief Evolve pseudospectral equations
          */
         virtual void evolveAfterPrognostic(const bool finishedStep);

         /**
          * @brief Add scalar equation to solver
          *
          * @param spEq Equation
          * @param it   Iteration index
          */
         void addEquation(Equations::SharedIScalarEquation spEq, const int it);

         /**
          * @brief Add scalar equation to solver
          *
          * @param spEq       Equation
          * @param eqId       Equation ID
          * @param it         Iteration index
          */
         void addEquation(Equations::SharedIScalarEquation spEq, const std::size_t eqId, const int it);

         /**
          * @brief Add vector equation to solver
          *
          * @param spEq Equation
          * @param it   Iteration index
          */
         void addEquation(Equations::SharedIVectorEquation spEq, const int it);

         /**
          * @brief Add vector equation to solver
          *
          * @param spEq       Equation
          * @param eqId       Equation ID
          * @param it         Iteration index
          */
         void addEquation(Equations::SharedIVectorEquation spEq, const std::size_t eqId, const int it);

         /**
          * @brief Add graph description to solver
          */
         void addGraph(const std::string& graphStr);

         /**
          * @brief Initialise the transforms
          *
          * @param descr Description of parallelisation
          */
         void initParallel(SharedResolution spRes, const Parallel::SplittingDescription& descr);

         /**
          * @brief Initialise the different components of the simulation
          *
          * @param tstep Configuration timestep information
          * @param spBcs Boundary condition information
          */
         void init(const Array& tstep, const SharedSimulationBoundary spBcs);

         /**
          * @brief Use state file time and timestep for diagnostics
          */
         void useStateTime(const MHDFloat time, const MHDFloat timestep);

         /**
          * @brief Scalar variables
          */
         const ScalarVariable_map& scalarVariables();

         /**
          * @brief Vector variables
          */
         const VectorVariable_map& vectorVariables();

         /**
          * @brief Imposed scalar variables
          */
         const ScalarVariable_map& imposedScalarVariables();

         /**
          * @brief Imposed vector variables
          */
         const VectorVariable_map& imposedVectorVariables();

         /**
          * @brief Transform coordinator
          */
         Transform::TransformCoordinatorType& transformCoordinator();

         /**
          * @brief Initialise the solvers (done just before preRun)
          */
         virtual void initSolvers();

         /**
          * @brief Last cleanup before run
          */
         virtual void cleanupForRun();

         /**
          * @brief Update physical space fields
          *
          * @param it Iteration index
          */
         void updatePhysical(const int it);

         /**
          * @brief Update spectral space fields
          *
          * @param it Iteration index
          */
         void updateSpectral(const int it);

         /**
          * @brief Update spectral space fields
          *
          * @param isTrivial        Update trivial equations?
          * @param isDiagnostic     Update diagnostic equations?
          * @param isPrognostic     Update prognostic equations?
          * @param isWrapper        Update wrapper equations?
          * @param it               Iteration index
          */
         void updateSpectral(const bool isTrivial, const bool isDiagnostic, const bool isPrognostic, const bool isWrapper, const int it);

         /**
          * @brief Compute the nonlinear terms
          *
          * @param it   Iteration index
          */
         void computeNonlinear(const int it);

         /**
          * @brief Explicit linear the trivial equations
          */
         void explicitTrivialEquations(const std::size_t opId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range);

         /**
          * @brief Explicit linear the trivial equations
          *
          * @param opId Operator ID
          * @param it   Iteration index
          */
         void explicitTrivialEquations(const std::size_t opId, const int it);

         /**
          * @brief Explicit linear for the diagnostic equations
          */
         void explicitDiagnosticEquations(const std::size_t opId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range);

         /**
          * @brief Explicit linear for the diagnostic equations
          *
          * @param opId Operator ID
          * @param it   Iteration index
          */
         void explicitDiagnosticEquations(const std::size_t opId, const int it);

         /**
          * @brief Solve the trivial equations
          */
         void solveTrivialEquations(const std::size_t timeId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range);

         /**
          * @brief Solve the trivial equations
          *
          * @param timeId  Timing ID
          * @param it      Iteration index
          */
         void solveTrivialEquations(const std::size_t timeId, const int it);

         /**
          * @brief Solve the diagnostic equations
          */
         void solveDiagnosticEquations(const std::size_t timeId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range);

         /**
          * @brief Solve the diagnostic equations
          *
          * @param timeId  Timing ID
          * @param it      Iteration index
          */
         void solveDiagnosticEquations(const std::size_t timeId, const int it);

         /**
          * @brief Explicit linear term for all equations
          *
          * @param it   Iteration index
          */
         void explicitEquations(const int it);

         /**
          * @brief Prepare time evolution
          *
          * @param schemeId ID of timestepping scheme
          */
         void prepareEvolution(const std::size_t schemeId);

         /**
          * @brief Pre solve equations for full initialisation
          */
         void preSolveEquations();

         /**
          * @brief Solve equations before prognostic equations
          *
          * @param it   Iteration index
          */
         void solveEquationsBefore(const int it);

         /**
          * @brief Solve equations after prognostic equations
          *
          * @param it   Iteration index
          */
         void solveEquationsAfter(const int it);

         /**
          * @brief Update equations
          *
          * @param it         Iteration index
          * @param isFinished full timestep computed?
          */
         void updateEquations(const int it, const bool isFinished) const;

         /**
          * @brief Finalize timestep iteration
          */
         void finalizeTimestep();

         /**
          * @brief Update the time stored in each equation
          */
         void updateEquationTime(const MHDFloat time, const bool finished);

         /**
          * @brief Explicit linear for the prognostic equations
          */
         void explicitPrognosticEquations(const std::size_t opId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range);

         /**
          * @brief Explicit linear for the prognostic equations
          *
          * @param it   Iteration index
          */
         void explicitPrognosticEquations(const std::size_t opId, const int it);

         /**
          * @brief Timestep the prognostic equations
          */
         void solvePrognosticEquations(ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range);

         /**
          * @brief Timestep the prognostic equations
          *
          * @param it   Iteration index
          */
         void solvePrognosticEquations(const int it);

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const;

         /**
          * @brief Profile storage used by pseudospectral coordinator
          */
         virtual void profileStorage() const;

         /**
          * @brief Write equation diagnostics
          *
          * @param isAsciiTime   ASCII writing time?
          * @param isHdf5Time    HDF5 writing time?
          */
         void writeDiagnostics(const bool isAsciiTime, const bool isHdf5Time) const;

      protected:
         /// Typedef for map of equations indexed by type ID and iteration index
         template <typename TEquation> using EquationTypeMapType = std::map<std::pair<std::size_t,int>, std::vector<std::shared_ptr<TEquation> > >;

         /// Typedef for map of equations indexed by iteration
         template <typename TEquation> using EquationIterMapType = std::map<int, std::vector<std::shared_ptr<TEquation> > >;

         /**
          * @brief Get scalar equation range
          *
          * @param eqId Equation ID
          * @param it   Iteration index
          */
         ScalarEquation_range scalarRange(const std::size_t eqId, const int it);

         /**
          * @brief Get vector equation range
          *
          * @param eqId Equation ID
          * @param it   Iteration index
          */
         VectorEquation_range vectorRange(const std::size_t eqId, const int it);

         /**
          * @brief Get scalar equation range
          *
          * @param eqId Equation ID
          * @param it   Iteration index
          */
         SharedIScalarEquation scalarEq(const std::size_t eqId, const int it, const int i = 0);

         /**
          * @brief Get vector equation range
          *
          * @param eqId Equation ID
          * @param it   Iteration index
          */
         SharedIVectorEquation vectorEq(const std::size_t eqId, const int it, const int i = 0);

         /**
          * @brief Initialise the imposed components of the simulation
          */
         void initImposed();

         /**
          * @brief Evolve pseudospectral equations
          *
          * @param it Current iteration
          */
         void evolveBefore(const int it);

         /**
          * @brief Evolve pseudospectral equations
          *
          * @param it Current iteration
          */
         void evolveEndIteration(const int it);

         /**
          * @brief Iteration indexes
          */
         std::set<int> mIt;

         /**
          * @brief Shared resolution
          */
         SharedResolution mspRes;

         /**
          * @brief Shared resolution
          */
         Equations::SharedCEquationParameters mspEqParams;

         /**
          * @brief Scalar equation ranges
          */
         EquationTypeMapType<IScalarEquation> mScalarEqMap;

         /**
          * @brief Vector equation ranges
          */
         EquationTypeMapType<IVectorEquation> mVectorEqMap;

         /**
          * @brief Trivial solver coordinator
          */
         SparseTrivialCoordinator mTrivialCoordinator;

         /**
          * @brief Linear solver coordinator
          */
         SparseLinearCoordinator mLinearCoordinator;

         /**
          * @brief Timestep coordinator
          */
         Timestep::Coordinator mTimestepCoordinator;

         /**
          * @brief Diagnostic coordinator
          */
         Diagnostics::Coordinator  mDiagnostics;

         /**
          * @brief Transform coordinator
          */
         Transform::TransformCoordinatorType mTransformCoordinator;

         /**
          * @brief Storage for scalar equations
          */
         EquationIterMapType<IScalarEquation> mScalarEquations;

         /**
          * @brief Storage for vector equations
          */
         EquationIterMapType<IVectorEquation> mVectorEquations;

         /**
          * @brief Map between name and pointer for the scalar variables
          */
         ScalarVariable_map  mScalarVariables;

         /**
          * @brief Map between name and pointer for the vector variables
          */
         VectorVariable_map  mVectorVariables;

         /**
          * @brief Transform coordinator
          */
         std::shared_ptr<Transform::TransformCoordinatorType> mspImposedTransformCoordinator;

         /**
          * @brief Map between name and pointer for the imposed scalar variables
          */
         ScalarVariable_map  mImposedScalarVariables;

         /**
          * @brief Map between name and pointer for the imposed vector variables
          */
         VectorVariable_map  mImposedVectorVariables;

         /**
          * @brief Map between name and physical kernels for each iteration
          */
         std::map<int, std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel> > mPhysicalKernels;

         /**
          * @brief Backward transform trees
          */
         std::map<int, std::vector<Transform::TransformTree> > mBwdTree;

         /**
          * @brief Forward transform trees
          */
         std::map<int, std::vector<Transform::TransformTree> > mFwdTree;

         /**
          * @brief Storage for a shared forward transform grouper
          */
         Transform::SharedIForwardGrouper   mspFwdGrouper;

         /**
          * @brief Storage for a shared backward transform grouper
          */
         Transform::SharedIBackwardGrouper   mspBwdGrouper;

         /**
          * @brief Storage for a shared forward transform grouper for imposed fields
          */
         Transform::SharedIForwardGrouper   mspImposedFwdGrouper;

         /**
          * @brief Storage for a shared backward transform grouper for imposed fields
          */
         Transform::SharedIBackwardGrouper   mspImposedBwdGrouper;


         std::unique_ptr<Graph::Jit<3>> mJitter;
         std::shared_ptr<QuICC::Memory::memory_resource> mMemRsr;

      private:
         /**
          * @brief If at least one equation is present in range
          */
         bool atLeastOne(ScalarEquation_range s, VectorEquation_range v) const;

         /**
          * @brief Iteration has prognostic equation
          */
         bool hasPrognostic(const int it) const;

         /**
          * @brief Add equation range to equation map
          */
         void addToMap(const std::size_t eqId, const int it, const ScalarEquation_range& r);

         /**
          * @brief Add equation range to equation map
          */
         void addToMap(const std::size_t eqId, const int it, const VectorEquation_range& r);

         /**
          * @brief Initialise the transform coordinator
          */
         void initTransformCoordinator();

         /**
          * @brief Initialise the equations (generate operators, etc)
          */
         void setupEquations();

         /**
          * @brief Sort equations and store information for timestep/solver/nothing ranges
          */
         void sortEquations();

         /**
          * @brief Setup the output files added by the model
          */
         void setupOutput();
   };

   inline const Coordinator::ScalarVariable_map& Coordinator::scalarVariables()
   {
      return this->mScalarVariables;
   }

   inline const Coordinator::VectorVariable_map& Coordinator::vectorVariables()
   {
      return this->mVectorVariables;
   }

   inline const Coordinator::ScalarVariable_map& Coordinator::imposedScalarVariables()
   {
      return this->mImposedScalarVariables;
   }

   inline const Coordinator::VectorVariable_map& Coordinator::imposedVectorVariables()
   {
      return this->mImposedVectorVariables;
   }

   inline Transform::TransformCoordinatorType& Coordinator::transformCoordinator()
   {
      return this->mTransformCoordinator;
   }

   /// Typedef for a shared pointer of a Coordinator
   typedef std::shared_ptr<Coordinator> SharedCoordinator;

}
}

#endif // QUICC_PSEUDOSPECTRAL_COORDINATOR_HPP
