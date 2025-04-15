/**
 * @file Interface.hpp
 * @brief Implementation of an interface to the ImEx predictor-corrector schemes
 */

#ifndef QUICC_TIMESTEP_PREDICTORCORRECTOR_INTERFACEVIEWS_HPP
#define QUICC_TIMESTEP_PREDICTORCORRECTOR_INTERFACEVIEWS_HPP

// System includes
//
#include <memory>
#include <type_traits>

// Project includes
//
#include "Memory/MemoryResource.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Pseudospectral/Coordinator.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/Timestep/Constants.hpp"
#include "QuICC/Timestep/IScheme.hpp"
#include "QuICC/Timestep/Interface.hpp"
#include "Timestep/PredictorCorrector/Views/TimestepperCoordinator.hpp"
#include "Timestep/PredictorCorrector/Views/ImExPCTimestepper.hpp"
#include "Timestep/PredictorCorrector/Views/details/TimesteppperTools.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "View/ViewDense.hpp"
#include "Memory/Memory.hpp"
#include "Memory/Cpu/NewDelete.hpp"

namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

/**
 * @brief Implementation of general timestepper structure
 */
template <typename TScheme> class InterfaceViews : public Timestep::Interface
{
public:
   /// Typedef for solver implementation
   template <typename T1, typename T2, typename T3>
   using SolverImplementationType =
      Views::ImExPCTimestepper<T1, T2, T3>;

   /// Typedef for parent coordinator
   typedef Views::TimestepperCoordinator<SolverImplementationType, base_t>
      SolverCoordinator;

   /**
    * @brief Constructor
    *
    * @param time    Initial time value
    * @param cfl     Initial CFL timestep
    * @param error   Max error allowed during timestep
    * @param scalEq  Shared scalar equations
    * @param vectEq  Shared vector equations
    */
   InterfaceViews(const MHDFloat time, const Matrix& cfl, const MHDFloat maxError,
      const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq,
      Pseudospectral::Coordinator& pseudo);

   /**
    * @brief Destructor
    */
   virtual ~InterfaceViews() = default;

   /**
    * @brief Timestep is finished?
    */
   bool finishedStep() const final;

   /**
    * @brief Set solve time
    */
   void setSolveTime(const std::size_t timeId) final;

   /**
    * @brief Update equation explicit linear input to solver
    *
    * @param scalEq Scalar equations
    * @param vectEq Vector equations
    */
   void getExplicitInput(const std::size_t opId,
      const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq,
      const ScalarVariable_map& scalVar,
      const VectorVariable_map& vectVar) final;

   /**
    * @brief Tune adaptive timestepper
    *
    * \mhdBug Not fully implemented
    */
   void tuneAdaptive(const MHDFloat time) final;

   /**
    * @brief Adapt the timestep used
    *
    * @param cfl     CFL conditions
    * @param scalEq  Shared scalar equations
    * @param vectEq  Shared vector equations
    */
   void adaptTimestep(const Matrix& cfl, const ScalarEquation_range& scalEq,
      const VectorEquation_range& vectEq) final;

   /**
    * @brief Compute (partial) forward step
    *
    * @param scalEq Shared scalar equations
    * @param vectEq Shared vector equations
    * @param scalVar Shared scalar variables
    * @param vectVar Shared vector variables
    */
   void stepForward(const ScalarEquation_range& scalEq,
      const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar,
      const VectorVariable_map& vectVar) final;

   /**
    * @brief Print timestepper information to stream
    *
    * @param stream  Output stream
    */
   void printInfo(std::ostream& stream) final;

protected:
   /**
    * @brief Initialize solution
    *
    * @param scalEq Scalar equations
    * @param vectEq Vector equations
    */
   void initSolution(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

   /**
    * @brief Update equation input to solver
    *
    * @param scalEq Scalar equations
    * @param vectEq Vector equations
    */
   void getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

   /**
    * @brief Transfer solution from solver
    *
    * @param scalEq Scalar equations
    * @param vectEq Vector equations
    */
   void transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

   /**
    * @brief Translate equations to timestepper info
    *
    * @param infos  Vector of information
    * @param scalEq Scalar equations
    * @param vectEq Vector equations
    */
   void translate(std::vector<Views::TimestepperInfo>& infos,
      const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

private:
   /**
    * @brief Interface to timestepping scheme
    */
   SharedIScheme mspScheme;

   /**
    * @brief Interface to timestepping scheme
    */
   SolverCoordinator mSolverCoord;

   /**
    * @brief
    */
   std::shared_ptr<Memory::memory_resource> _mem;
};

/**
 * @brief Compute flat timestepper index
 */
std::size_t stepperIndex(const std::size_t solverIndex, const std::size_t sysIndex);


/**
 * @brief Wrapper to build timestepping matrices
 */
void buildTimestepMatrixWrapper(std::map<std::size_t, DecoupledZSparse>& ops, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

template <typename TScheme>
InterfaceViews<TScheme>::InterfaceViews(const MHDFloat time, const Matrix& cfl,
   const MHDFloat maxError, const ScalarEquation_range& scalEq,
   const VectorEquation_range& vectEq, Pseudospectral::Coordinator& pseudo) :
    Timestep::Interface(time, cfl, maxError, scalEq, vectEq, pseudo)
{
   _mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();

   std::cerr << "##############################################" << std::endl;
   std::cerr << "##############################################" << std::endl;
   std::cerr << "##############################################" << std::endl;
   std::cerr << "INITIALIZING NEW TIMESTEPPER INFRASTRUCTURE" << std::endl;
   std::cerr << "##############################################" << std::endl;
   std::cerr << "##############################################" << std::endl;
   std::cerr << "##############################################" << std::endl;

   // Create Timestepper scheme
   std::shared_ptr<TScheme> spScheme = std::make_shared<TScheme>();

   // Use embedded scheme to compute error
   if (maxError > 0.0)
   {
      spScheme->enableEmbedded();
      this->mMaxError = maxError;
   }
   this->mspScheme = spScheme;

   std::vector<Views::TimestepperInfo> infos;
   this->translate(infos, scalEq, vectEq);

   this->mSolverCoord.init(this->timestep(), infos, spScheme);

   this->initSolution(scalEq, vectEq);
}

// \todo convert to free function
template <typename TScheme>
void InterfaceViews<TScheme>::translate(std::vector<Views::TimestepperInfo>& infos, const ScalarEquation_range& scalEq,
   const VectorEquation_range& vectEq)
{
   auto addInfo = [](auto& infos, auto&& eq_range)
   {
      SpectralFieldId myId;

      // Loop over range
      for(auto& eqIt: make_range(eq_range))
      {
         for(auto& compIt: make_range(eqIt->spectralRange()))
         {
            // Get field identity
            myId = std::make_pair(eqIt->name(), compIt);

            const auto& cinfo = eqIt->couplingInfo(myId.second);
            typename std::remove_reference_t<decltype(infos)>::value_type info;
            info.isComplex = cinfo.isComplex();
            info.fieldIndex = cinfo.fieldIndex();
            info.timeId = eqIt->solveTiming();
            for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
            {
               info.solverIndex = stepperIndex(cinfo.solverIndex(), i);
               info.rows = cinfo.systemN(i);
               info.cols = cinfo.rhsCols(i);
               info.blockN = cinfo.galerkinN(i);
               // Set operators
               buildTimestepMatrixWrapper(info.ops, eqIt, myId.second, i);

               infos.push_back(info);
            }
         }
      }
   };

   addInfo(infos, scalEq);
   addInfo(infos, vectEq);
}

inline std::size_t stepperIndex(const std::size_t solverIndex, const std::size_t sysIndex)
{
   assert(sysIndex < 10000);

   return 10000*solverIndex + sysIndex;
}

template <typename TScheme>
void InterfaceViews<TScheme>::tuneAdaptive(const MHDFloat time)
{
   this->mStepTime = time;
}

template <typename TScheme> bool InterfaceViews<TScheme>::finishedStep() const
{
   return this->mSolverCoord.finishedStep();
}

template <typename TScheme>
void InterfaceViews<TScheme>::setSolveTime(const std::size_t timeId)
{
   this->mSolverCoord.setSolveTime(timeId);
}

template <typename TScheme>
void InterfaceViews<TScheme>::initSolution(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
{
   auto processSolution = [this](const auto currentTime, auto&& eq_range)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      for(auto& eqIt: make_range(eq_range))
      {
         if(eqIt->solveTiming() == currentTime)
         {
            for(auto& compIt: make_range(eqIt->spectralRange()))
            {
               // Get field identity
               myId = std::make_pair(eqIt->name(), compIt);

               const auto& cinfo = eqIt->couplingInfo(myId.second);
               Views::TimestepperInfo info;
               info.isComplex = cinfo.isComplex();
               info.fieldIndex = cinfo.fieldIndex();
               info.timeId = eqIt->solveTiming();

               for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
               {
                  DecoupledZMatrix tmp(cinfo.galerkinN(i), cinfo.rhsCols(i));
                  tmp.setZero();

                  std::visit(
                        [&](auto&& p)
                        {
                        Equations::copyUnknown(*eqIt, p->dom(0).perturbation(), myId.second, tmp, i, 0, true, true);
                        }, eqIt->spUnknown());

                  std::uint32_t mem_rows = static_cast<std::uint32_t>(cinfo.galerkinN(i));
                  std::uint32_t mem_cols = static_cast<std::uint32_t>(cinfo.rhsCols(i));
                  Memory::MemBlock<MHDComplex> data(mem_rows*mem_cols, this->_mem.get());
                  using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
                  std::array<std::uint32_t, 2> dimensions {mem_rows, mem_cols};
                  View::View<MHDComplex, View::Attributes<dense2D>> tmpView(data, dimensions); 
                  Views::details::computeSet(tmpView, tmp);

                  this->mSolverCoord.updateSolution(info, tmpView);
               }
            }
         }
      }
   };

   processSolution(this->mSolverCoord.solveTime(), scalEq);
   processSolution(this->mSolverCoord.solveTime(), vectEq);
}

template <typename TScheme>
void InterfaceViews<TScheme>::getExplicitInput(const std::size_t opId,
   const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq,
   const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
{
   Profiler::RegionFixture<2> fix("Timestep-explicitInput");

   auto processInput = [this](const auto currentTime, auto&& eq_range, const std::size_t opId,
   const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      for(auto& eqIt: make_range(eq_range))
      {
         if(eqIt->solveTiming() == currentTime)
         {
            for(auto& compIt: make_range(eqIt->spectralRange()))
            {
               // Get field identity
               myId = std::make_pair(eqIt->name(), compIt);

               // Apply constraint on solution
               eqIt->applyConstraint(myId.second, SolveTiming::Before::id());

               const auto& cinfo = eqIt->couplingInfo(myId.second);
               Views::TimestepperInfo info;
               info.isComplex = cinfo.isComplex();
               info.fieldIndex = cinfo.fieldIndex();
               info.timeId = eqIt->solveTiming();

               // Get timestep input
               for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
               {
                  info.solverIndex = stepperIndex(cinfo.solverIndex(), i);
                  info.rows = cinfo.systemN(i);
                  info.cols = cinfo.rhsCols(i);
                  info.blockN = cinfo.galerkinN(i);

                  // Copy field values into solver input
                  DecoupledZMatrix tmp(cinfo.galerkinN(i), cinfo.rhsCols(i));
                  tmp.setZero();

                  // Build range of operator
                  auto r = make_range(cinfo.explicitRange(opId));

#ifdef QUICC_DEBUG
                  if(r.size() == 0)
                  {
                     DebuggerMacro_msg("(Nothing)", 7);
                  }
#endif // QUICC_DEBUG

                  // Loop over explicit fields
                  for(auto& fIt: r)
                  {
                     DebuggerMacro_msg("Add " + ModelOperator::Coordinator::tag(opId) + " term from " + PhysicalNames::Coordinator::tag(fIt.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(fIt.second)) + ")", 7);

                     // Get explicit input
                     if(fIt.second == FieldComponents::Spectral::SCALAR)
                     {
                        std::visit(
                              [&](auto&& p)
                              {
                              Equations::addExplicitTerm(*eqIt, opId, myId.second, tmp, 0, fIt, p->dom(0).perturbation(), i);
                              }, scalVar.find(fIt.first)->second);
                     } else
                     {
                        std::visit(
                              [&](auto&& p)
                              {
                              Equations::addExplicitTerm(*eqIt, opId, myId.second, tmp, 0, fIt, p->dom(0).perturbation().comp(fIt.second), i);
                              }, vectVar.find(fIt.first)->second);
                     }
                  }

                  std::uint32_t mem_rows = static_cast<std::uint32_t>(cinfo.galerkinN(i));
                  std::uint32_t mem_cols = static_cast<std::uint32_t>(cinfo.rhsCols(i));
                  Memory::MemBlock<MHDComplex> data(mem_rows*mem_cols, this->_mem.get());
                  using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
                  std::array<std::uint32_t, 2> dimensions {mem_rows, mem_cols};
                  View::View<MHDComplex, View::Attributes<dense2D>> tmpView(data, dimensions); 
                  Views::details::computeSet(tmpView, tmp);

                  this->mSolverCoord.updateRhs(info, tmpView);
               }
            }
         }
      }
   };

   processInput(this->mSolverCoord.solveTime(), scalEq, opId, scalVar, vectVar);
   processInput(this->mSolverCoord.solveTime(), vectEq, opId, scalVar, vectVar);
}

template <typename TScheme>
void InterfaceViews<TScheme>::getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
{
   Profiler::RegionFixture<2> fix("Timestep-input");

   auto processInput = [this](const auto currentTime, auto&& eq_range)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      for(auto& eqIt: make_range(eq_range))
      {
         if(eqIt->solveTiming() == currentTime)
         {
            for(auto& compIt: make_range(eqIt->spectralRange()))
            {
               // Get field identity
               myId = std::make_pair(eqIt->name(), compIt);

               // Apply constraint on solution
               eqIt->applyConstraint(myId.second, SolveTiming::Before::id());

               const auto& cinfo = eqIt->couplingInfo(myId.second);
               Views::TimestepperInfo info;
               info.isComplex = cinfo.isComplex();
               info.fieldIndex = cinfo.fieldIndex();
               info.timeId = eqIt->solveTiming();

               // Get timestep input
               for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
               {
                  info.solverIndex = stepperIndex(cinfo.solverIndex(), i);
                  info.rows = cinfo.systemN(i);
                  info.cols = cinfo.rhsCols(i);
                  info.blockN = cinfo.galerkinN(i);

                  // Copy field values into solver input
                  DecoupledZMatrix tmp(cinfo.galerkinN(i), cinfo.rhsCols(i));
                  tmp.setZero();
                  Equations::copyNonlinear(*eqIt, myId.second, tmp, i, 0);
                  std::cerr << "COPY NONLINEAR DATA: Re" << std::endl;
                  std::cerr << tmp.real() << std::endl;
                  std::cerr << "COPY NONLINEAR DATA: Im" << std::endl;
                  std::cerr << tmp.imag() << std::endl;

                  // Add source term
                  std::visit(
                        [&](auto&& p)
                        {
                        Equations::addSource(*eqIt, p->dom(0).perturbation(), myId.second, tmp, i, 0);
                        }, eqIt->spUnknown());

                  std::uint32_t mem_rows = static_cast<std::uint32_t>(cinfo.galerkinN(i));
                  std::uint32_t mem_cols = static_cast<std::uint32_t>(cinfo.rhsCols(i));
                  Memory::MemBlock<MHDComplex> data(mem_rows*mem_cols, this->_mem.get());
                  using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
                  std::array<std::uint32_t, 2> dimensions {mem_rows, mem_cols};
                  View::View<MHDComplex, View::Attributes<dense2D>> tmpView(data, dimensions); 
                  Views::details::computeSet(tmpView, tmp);

                  this->mSolverCoord.updateRhs(info, tmpView);

                  // If required set inhomogenous boundary condition value
                  if(cinfo.hasBoundaryValue())
                  {
                     // Set boundary value
                     std::visit(
                           [&](auto&& p)
                           {
                              throw std::logic_error("NOT YET IMPLEMENTED");
                             // Equations::setBoundaryValue(*spEq, p->dom(0).perturbation(), id.second, (*solveIt)->rInhomogeneous(i), i, (*solveIt)->startRow(id,i));
                           }, eqIt->spUnknown());
      //               this->mSolverCoord.updateInhomogeneous(info);
                  }
               }
            }
         }
      }
   };

   processInput(this->mSolverCoord.solveTime(), scalEq);
   processInput(this->mSolverCoord.solveTime(), vectEq);
}

template <typename TScheme>
void InterfaceViews<TScheme>::transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
{
   Profiler::RegionFixture<2> fix("Timestep-output");

   auto processOutput = [this](const auto currentTime, auto&& eq_range)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      for(auto& eqIt: make_range(eq_range))
      {
         if(eqIt->solveTiming() == currentTime)
         {
            for(auto& compIt: make_range(eqIt->spectralRange()))
            {
               // Get field identity
               myId = std::make_pair(eqIt->name(), compIt);

               const auto& cinfo = eqIt->couplingInfo(myId.second);
               Views::TimestepperInfo info;
               info.isComplex = cinfo.isComplex();
               info.fieldIndex = cinfo.fieldIndex();
               info.timeId = eqIt->solveTiming();

               for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
               {
                  info.solverIndex = stepperIndex(cinfo.solverIndex(), i);
                  info.rows = cinfo.systemN(i);
                  info.cols = cinfo.rhsCols(i);
                  info.blockN = cinfo.galerkinN(i);

                  std::uint32_t mem_rows = static_cast<std::uint32_t>(cinfo.galerkinN(i));
                  std::uint32_t mem_cols = static_cast<std::uint32_t>(cinfo.rhsCols(i));
                  Memory::MemBlock<MHDComplex> data(mem_rows*mem_cols, this->_mem.get());
                  using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
                  std::array<std::uint32_t, 2> dimensions {mem_rows, mem_cols};
                  View::View<MHDComplex, View::Attributes<dense2D>> tmpView(data, dimensions); 
                  this->mSolverCoord.getSolution(tmpView, info);

                  DecoupledZMatrix tmp(cinfo.galerkinN(i), cinfo.rhsCols(i));
                  Views::details::computeSet(tmp, tmpView);

                  eqIt->storeSolution(myId.second, tmp, i, 0);
               }

               // Apply constraint on solution
               auto changedSolution = eqIt->applyConstraint(myId.second, SolveTiming::After::id());

               // Update solver solution if constraint modified it
               if(changedSolution)
               {
                  for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
                  {
                     info.solverIndex = stepperIndex(cinfo.solverIndex(), i);
                     info.rows = cinfo.systemN(i);
                     info.cols = cinfo.rhsCols(i);
                     info.blockN = cinfo.galerkinN(i);

                     DecoupledZMatrix tmp(cinfo.galerkinN(i), cinfo.rhsCols(i));
                     tmp.setZero();

                     std::visit(
                           [&](auto&& p)
                           {
                           Equations::copyUnknown(*eqIt, p->dom(0).perturbation(), myId.second, tmp, i, 0, true, true);
                           }, eqIt->spUnknown());

                     std::uint32_t mem_rows = static_cast<std::uint32_t>(cinfo.galerkinN(i));
                     std::uint32_t mem_cols = static_cast<std::uint32_t>(cinfo.rhsCols(i));
                     Memory::MemBlock<MHDComplex> data(mem_rows*mem_cols, this->_mem.get());
                     using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
                     std::array<std::uint32_t, 2> dimensions {mem_rows, mem_cols};
                     View::View<MHDComplex, View::Attributes<dense2D>> tmpView(data, dimensions); 
                     Views::details::computeSet(tmpView, tmp);

                     this->mSolverCoord.updateSolution(info, tmpView);
                  }
               }
            }
         }
      }
   };

   processOutput(this->mSolverCoord.solveTime(), scalEq);
   processOutput(this->mSolverCoord.solveTime(), vectEq);
}

template <typename TScheme>
void InterfaceViews<TScheme>::adaptTimestep(const Matrix& cfl,
   const ScalarEquation_range&, const VectorEquation_range&)
{
   // Store old timestep
   this->mOldDt = this->timestep();

   // Update CFL information
   this->mDt.block(0, 1, this->mDt.rows(), cfl.cols() - 1) =
      cfl.rightCols(cfl.cols() - 1);

   // New computed CFL
   MHDFloat compCfl = cfl(0, 0);

   // Check if CFL allows for a larger timestep
   MHDFloat newCflDt = 0.0;
   if (compCfl > this->mcUpWindow * this->timestep())
   {
      if (this->mCnstSteps >= this->mcMinCnst)
      {
         // Set new timestep
         newCflDt = std::min(compCfl, this->mcMaxJump * this->timestep());
      }
      else
      {
         // Reuse same timestep
         newCflDt = this->timestep();
      }

      // Check if CFL is below minimal timestep or downard jump is large
   }
   else if (compCfl < this->mcMinDt ||
            compCfl < this->timestep() / this->mcMaxJump)
   {
      // Signal simulation abort
      newCflDt = -compCfl;

      // Check if CFL requires a lower timestep
   }
   else if (compCfl < this->timestep() * (2.0 - this->mcUpWindow))
   {
      // Set new timestep
      newCflDt = compCfl;
   }
   else
   {
      newCflDt = this->timestep();
   }

   // Get timestepper error (if applicable)
   MHDFloat error = this->mSolverCoord.error();

// Gather error across processes
#ifdef QUICC_MPI
   if (error > 0.0)
   {
      MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX,
         MPI_COMM_WORLD);
   }
#endif // QUICC_MPI

   // No error control and no CFL condition
   MHDFloat newErrorDt = 0.0;

   // Use what ever condition is used by CFL
   if (error < 0)
   {
      newErrorDt = -1.0;

      // Error is too large, reduce timestep
   }
   else if (error > this->mMaxError)
   {
      newErrorDt =
         this->timestep() *
         std::pow(this->mMaxError / error, 1. / this->mspScheme->order()) /
         this->mcUpWindow;

      // Error is small, increase timestep
   }
   else if (error < this->mMaxError / (this->mcMaxJump * 0.9) &&
            this->mCnstSteps >= this->mcMinCnst)
   {
      newErrorDt =
         std::min(this->timestep() * std::pow(this->mMaxError / error,
                                        1. / this->mspScheme->order()),
            this->timestep() * this->mcMaxJump);

      // Timestep should not be increased
   }
   else
   {
      newErrorDt = this->timestep();
   }

   // Update error details
   if (this->mMaxError > 0.0)
   {
      this->mDt(0, this->mDt.cols() - 1) = newErrorDt;
      this->mDt(1, this->mDt.cols() - 1) = error;
   }

   // CFL condition requested abort!
   if (newCflDt < 0.0)
   {
      this->mDt(0, 0) = newCflDt;
      this->mDt(1, 0) = cfl(1, 0);

      // Get minimum between both conditions
   }
   else if (newCflDt > 0.0 && newErrorDt > 0.0)
   {
      if (newCflDt < newErrorDt)
      {
         this->mDt(0, 0) = newCflDt;
         this->mDt(1, 0) = cfl(1, 0);
      }
      else
      {
         this->mDt(0, 0) = newErrorDt;
         this->mDt(1, 0) = ERROR_LOCATION;
      }

      // Use CFL condition
   }
   else if (newCflDt > 0.0)
   {
      if (this->timestep() != newCflDt)
      {
         this->mDt(0, 0) = newCflDt;
         this->mDt(1, 0) = cfl(1, 0);
      }

      // Use error condition
   }
   else if (newErrorDt > 0.0)
   {
      this->mDt(0, 0) = newErrorDt;
      this->mDt(1, 0) = ERROR_LOCATION;
   }

   //
   // Update the timestep matrices if necessary
   //
   if (this->timestep() != this->mOldDt && this->timestep() > 0.0)
   {
      DebuggerMacro_showValue(
         "Updating timestep and matrices with new Dt = ", 0, this->timestep());

      this->mSolverCoord.updateTimestep(this->timestep());

      // Update the time dependence in matrices
      DebuggerMacro_start("Update matrices", 0);
      this->mSolverCoord.updateMatrices();
      DebuggerMacro_stop("Update matrices t = ", 0);
   }
   else
   {
      this->mCnstSteps += 1.0;
   }

   // Update CFL writer
   this->mspIo->setSimTime(this->mTime, this->mDt, this->mCnstSteps);
   this->mspIo->write();

   if (this->timestep() != this->mOldDt && this->timestep() > 0.0)
   {
      this->mCnstSteps = 0.0;
   }
}

template <typename TScheme>
void InterfaceViews<TScheme>::stepForward(const ScalarEquation_range& scalEq,
   const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar,
   const VectorVariable_map& vectVar)
{
   bool isIntegrating = true;
   while (isIntegrating)
   {
      DebuggerMacro_msg("Time integration sub-step", 2);

      this->mpPseudo->evolveUntilPrognostic(this->finishedStep());
      this->setSolveTime(SolveTiming::Prognostic::id());

      // Update the equation input to the timestepper
      this->getInput(scalEq, vectEq);

      Profiler::RegionStart<2>("Timestep-solve");
      // Solve all the linear systems
      this->mSolverCoord.solveSystems();
      Profiler::RegionStop<2>("Timestep-solve");

      // Transfer timestep output back to equations
      this->transferOutput(scalEq, vectEq);

      // Clear the solver RHS
      this->mSolverCoord.clearSolvers();

      // Update current time
      this->mTime =
         this->mRefTime + this->mSolverCoord.stepFraction() * this->timestep();

      this->mpPseudo->evolveAfterPrognostic(this->finishedStep());

      isIntegrating = !this->finishedStep();
   }
}

template <typename TScheme>
void InterfaceViews<TScheme>::printInfo(std::ostream& stream)
{
   // Create nice looking ouput header
   Tools::Formatter::printNewline(stream);
   Tools::Formatter::printLine(stream, '-');
   Tools::Formatter::printCentered(stream, "Timestepper information", '*');
   Tools::Formatter::printLine(stream, '-');

   std::stringstream oss;
   int base = 20;

   // Timestep scheme
   oss << "Timestepper: " << this->mspScheme->name() << " ("
       << this->mspScheme->order() << ")";
   Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
   oss.str("");

   // General linear solver
   oss << "General solver: ";
#if defined QUICC_SPLINALG_MUMPS
   oss << "MUMPS";
#elif defined QUICC_SPLINALG_UMFPACK
   oss << "UmfPack";
#elif defined QUICC_SPLINALG_SPARSELU
   oss << "SparseLU";
#else
   oss << "(unknown)";
#endif // defined QUICC_SPLINALG_MUMPS

   Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
   oss.str("");

   // Triangular linear solver
   oss << "Triangular solver: ";
#if defined QUICC_SPTRILINALG_SPARSELU
   oss << "SparseLU";
#elif defined QUICC_SPTRILINALG_MUMPS
   oss << "MUMPS";
#elif defined QUICC_SPTRILINALG_UMFPACK
   oss << "UmfPack";
#else
   oss << "(unknown)";
#endif // defined QUICC_SPTRILINALG_SPARSELU

   Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
   oss.str("");

   // SPD linear solver
   oss << "SPD solver: ";
#if defined QUICC_SPSPDLINALG_SIMPLICIALLDLT
   oss << "SimplicialLDLT";
#elif defined QUICC_SPSPDLINALG_SIMPLICIALLLT
   oss << "SimplicialLLT";
#elif defined QUICC_SPSPDLINALG_MUMPS
   oss << "MUMPS";
#elif defined QUICC_SPSPDLINALG_UMFPACK
   oss << "UmfPack";
#elif defined QUICC_SPSPDLINALG_SPARSELU
   oss << "SparseLU";
#else
   oss << "(unknown)";
#endif // defined QUICC_SPSPDLINALG_SIMPLICIALLDLT

   Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
   oss.str("");

   Tools::Formatter::printLine(stream, '*');
   Tools::Formatter::printNewline(stream);
}

// \todo  move to separate file
inline void buildTimestepMatrixWrapper(std::map<std::size_t, DecoupledZSparse>& ops, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp,
   const int idx)
{
   bool isSplit = spEq->couplingInfo(comp).isSplitEquation();

   // Compute model's linear operator (without Tau lines)
   ops.insert(
      std::make_pair(ModelOperator::ImplicitLinear::id(), DecoupledZSparse()));
   spEq->buildModelMatrix(ops.find(ModelOperator::ImplicitLinear::id())->second,
      ModelOperator::ImplicitLinear::id(), comp, idx,
      ModelOperatorBoundary::SolverNoTau::id());
   // Compute model's time operator (without Tau lines)
   ops.insert(std::make_pair(ModelOperator::Time::id(), DecoupledZSparse()));
   spEq->buildModelMatrix(ops.find(ModelOperator::Time::id())->second,
      ModelOperator::Time::id(), comp, idx,
      ModelOperatorBoundary::SolverNoTau::id());
   // Compute model's tau line boundary operator
   ops.insert(
      std::make_pair(ModelOperator::Boundary::id(), DecoupledZSparse()));
   spEq->buildModelMatrix(ops.find(ModelOperator::Boundary::id())->second,
      ModelOperator::Boundary::id(), comp, idx,
      ModelOperatorBoundary::SolverHasBc::id());

   // If equation was split into two lower order systems
   if (isSplit)
   {
      // Compute model's split linear operator (without Tau lines)
      auto id = ModelOperator::SplitImplicitLinear::id();
      ops.insert(std::make_pair(id, DecoupledZSparse()));
      spEq->buildModelMatrix(ops.find(id)->second, id, comp, idx,
         ModelOperatorBoundary::SolverNoTau::id());

      // Compute model's tau line boundary operator for split operator
      id = ModelOperator::SplitBoundary::id();
      ops.insert(std::make_pair(id, DecoupledZSparse()));
      spEq->buildModelMatrix(ops.find(id)->second, id, comp, idx,
         ModelOperatorBoundary::SolverHasBc::id());

      // Compute model's tau line boundary value for split operator
      id = ModelOperator::SplitBoundaryValue::id();
      ops.insert(std::make_pair(id, DecoupledZSparse()));
      spEq->buildModelMatrix(ops.find(id)->second, id, comp, idx,
         ModelOperatorBoundary::SolverNoTau::id());
   }
}

} // namespace PredictorCorrector
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_PREDICTORCORRECTOR_INTERFACEVIEWS_HPP
