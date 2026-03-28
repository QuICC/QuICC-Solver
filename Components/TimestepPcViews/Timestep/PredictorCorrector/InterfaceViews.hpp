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
#include "QuICC/Equations/CopyUnknown.hpp"
#include "QuICC/ModelOperator/QuasiInverse.hpp"
#include "QuICC/ModelOperator/SplitQuasiInverse.hpp"
#include "QuICC/Equations/CorrectSolution.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Equations/AddSource.hpp"
#include "QuICC/Equations/SetBoundaryValue.hpp"
#include "QuICC/Equations/CopyNonlinear.hpp"
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

#define QUICC_DETAIL_PROF_LVL 3

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
    */
   void adaptTimestep(const Matrix& cfl) final;

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
   using Timestep::Interface::printInfo;

   /**
    * @brief Timestep is finished?
    */
   bool finishedStep() const;

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

Views::TimestepperInfo createInfo(const Equations::CouplingInformation& cinfo);

inline Views::TimestepperInfo createInfo(const Equations::CouplingInformation& cinfo, const std::size_t idx)
{
   Views::TimestepperInfo info;
   info.isComplex = cinfo.isComplex();
   info.fieldIndex = cinfo.fieldIndex();
   info.solverIndex = stepperIndex(cinfo.solverIndex(), idx);
   info.rows = cinfo.systemN(idx);
   info.cols = cinfo.rhsCols(idx);
   info.blockN = cinfo.galerkinN(idx);

   return info;
}

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
            if(eqIt->solveTiming() == SolveTiming::Prognostic::id())
            {
               DebuggerMacro_msg("Creating timesteppers for " + PhysicalNames::Coordinator::tag(eqIt->name()) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(compIt)) + ")", 2);

               for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
               {
                  auto info = createInfo(cinfo, i);

                  // Set operators
                  buildTimestepMatrixWrapper(info.ops, eqIt, myId.second, i);

                  infos.push_back(info);
               }
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
void InterfaceViews<TScheme>::initSolution(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
{
   auto processSolution = [this](auto&& eq_range)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      for(auto& eqIt: make_range(eq_range))
      {
            for(auto& compIt: make_range(eqIt->spectralRange()))
            {
               // Get field identity
               myId = std::make_pair(eqIt->name(), compIt);

               const auto& cinfo = eqIt->couplingInfo(myId.second);

               if(eqIt->solveTiming() == SolveTiming::Prognostic::id())
               {
                  // Allocate temporary storage
                  std::uint32_t mem_size = 0;
                  for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
                  {
                     mem_size = std::max(mem_size, static_cast<std::uint32_t>(cinfo.galerkinN(i))*static_cast<std::uint32_t>(cinfo.rhsCols(i)));
                  }
                  Memory::MemBlock<MHDComplex> data(mem_size, this->_mem.get());

                  for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
                  {
                     auto info = createInfo(cinfo, i);

                     std::uint32_t mem_rows = static_cast<std::uint32_t>(cinfo.galerkinN(i));
                     std::uint32_t mem_cols = static_cast<std::uint32_t>(cinfo.rhsCols(i));
                     using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
                     std::array<std::uint32_t, 2> dimensions {mem_rows, mem_cols};
                     View::View<MHDComplex, View::Attributes<dense2D>> tmpView(data, dimensions);

                     if(cinfo.isGalerkin())
                     {
                        std::visit(
                              [&](auto&& p)
                              {
                                 Equations::solveStencilUnknown(eqIt->res(), eqIt->couplingInfo(myId.second), myId.first, p->dom(0).perturbation(), myId.second, tmpView, i, 0, eqIt->backend(), eqIt->bcIds().map(), eqIt->eqParams().map());
                              }, eqIt->spUnknown());
                     }
                     else
                     {
                        std::visit(
                              [&](auto&& p)
                              {
                                 Equations::copyUnknown(eqIt->res(), eqIt->couplingInfo(myId.second), p->dom(0).perturbation(), myId.second, tmpView, i, 0, true, true, true);
                              }, eqIt->spUnknown());
                     }

                     this->mSolverCoord.updateSolution(info, tmpView);
                  }
               }
         }
      }
   };

   processSolution(scalEq);
   processSolution(vectEq);
}

template <typename TScheme>
void InterfaceViews<TScheme>::getExplicitInput(const std::size_t opId,
   const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq,
   const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
{
   Profiler::RegionFixture<2> fix("Timestep-explicitInput");

   auto processInput = [this](auto&& eq_range, const std::size_t opId,
   const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      for(auto& eqIt: make_range(eq_range))
      {
         for(auto& compIt: make_range(eqIt->spectralRange()))
         {
            // Get field identity
            myId = std::make_pair(eqIt->name(), compIt);

            const auto& cinfo = eqIt->couplingInfo(myId.second);

            if(eqIt->solveTiming() == SolveTiming::Prognostic::id())
            {
               DebuggerMacro_msg("Get explicit timestepper input for " + PhysicalNames::Coordinator::tag(myId.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(myId.second)) + ")", 6);

               // Build range of operator
               auto r = make_range(cinfo.explicitRange(opId));

#ifdef QUICC_DEBUG
               if(r.size() == 0)
               {
                  DebuggerMacro_msg("(Nothing)", 7);
               }
#endif // QUICC_DEBUG

               if(r.size() > 0)
               {
                  // Allocate temporary storage
                  std::uint32_t mem_size = 0;
                  for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
                  {
                     mem_size = std::max(mem_size, static_cast<std::uint32_t>(cinfo.galerkinN(i))*static_cast<std::uint32_t>(cinfo.rhsCols(i)));
                  }
                  Memory::MemBlock<MHDComplex> data(mem_size, this->_mem.get());

                  // Get timestep input
                  for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
                  {
                     auto info = createInfo(cinfo, i);

                     // Copy field values into timestepper input
                     DecoupledZMatrix tmp(cinfo.galerkinN(i), cinfo.rhsCols(i));
                     tmp.setZero();

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
                     using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
                     std::array<std::uint32_t, 2> dimensions {mem_rows, mem_cols};
                     View::View<MHDComplex, View::Attributes<dense2D>> tmpView(data, dimensions);
                     Views::details::computeSet(tmpView, tmp, 0);

                     this->mSolverCoord.updateRhs(info, tmpView);
                  }
               }
            }
         }
      }
   };

   processInput(scalEq, opId, scalVar, vectVar);
   processInput(vectEq, opId, scalVar, vectVar);
}

template <typename TScheme>
void InterfaceViews<TScheme>::getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
{
   Profiler::RegionFixture<2> fix("Timestep-input");

   auto processInput = [this](auto&& eq_range)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      for(auto& eqIt: make_range(eq_range))
      {
         for(auto& compIt: make_range(eqIt->spectralRange()))
         {
            // Get field identity
            myId = std::make_pair(eqIt->name(), compIt);

            // Apply constraint on solution
            eqIt->applyConstraint(myId.second, SolveTiming::Before::id());

            const auto& cinfo = eqIt->couplingInfo(myId.second);

            if(eqIt->solveTiming() == SolveTiming::Prognostic::id())
            {
   Profiler::RegionStart<QUICC_DETAIL_PROF_LVL>("Timestep-input:allocateTemp");
               // Allocate temporary storage
               std::uint32_t mem_size = 0;
               for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
               {
                  mem_size = std::max(mem_size, static_cast<std::uint32_t>(cinfo.galerkinN(i))*static_cast<std::uint32_t>(cinfo.rhsCols(i)));
               }
               Memory::MemBlock<MHDComplex> data(mem_size, this->_mem.get());
   Profiler::RegionStop<QUICC_DETAIL_PROF_LVL>("Timestep-input:allocateTemp");

               DebuggerMacro_msg("Get timestepper input for " + PhysicalNames::Coordinator::tag(myId.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(myId.second)) + ")", 6);

               // Get timestep input
               for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
               {
   Profiler::RegionStart<QUICC_DETAIL_PROF_LVL>("Timestep-input:createInfo");
                  auto info = createInfo(cinfo, i);
   Profiler::RegionStop<QUICC_DETAIL_PROF_LVL>("Timestep-input:createInfo");

   Profiler::RegionStart<QUICC_DETAIL_PROF_LVL>("Timestep-input:setupView");
                  std::uint32_t mem_rows = static_cast<std::uint32_t>(cinfo.galerkinN(i));
                  std::uint32_t mem_cols = static_cast<std::uint32_t>(cinfo.rhsCols(i));
                  using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
                  std::array<std::uint32_t, 2> dimensions {mem_rows, mem_cols};
                  View::View<MHDComplex, View::Attributes<dense2D>> tmpView(data, dimensions);
   Profiler::RegionStop<QUICC_DETAIL_PROF_LVL>("Timestep-input:setupView");

   Profiler::RegionStart<QUICC_DETAIL_PROF_LVL>("Timestep-input:copyNonlinear");
                  // Copy field values into timestepper input
                  if(cinfo.hasNonlinear())
                  {
                     std::visit(
                           [&](auto&& p)
                           {
                           Equations::copyUnknown(eqIt->res(), eqIt->couplingInfo(myId.second), p->dom(0).perturbation(), myId.second, tmpView, i, 0, true, true, false);
                           }, eqIt->spUnknown());
                  }
   Profiler::RegionStop<QUICC_DETAIL_PROF_LVL>("Timestep-input:copyNonlinear");

                  // Add source term
                  if(cinfo.hasSource())
                  {
                     std::visit(
                           [&](auto&& p)
                           {
                           Equations::addSource(eqIt->res(), cinfo, eqIt->sourceKernel(myId.second), p->dom(0).perturbation(), myId.second, tmpView, i, 0);
                           }, eqIt->spUnknown());
                  }

   Profiler::RegionStart<QUICC_DETAIL_PROF_LVL>("Timestep-input:updateRhs");
                  this->mSolverCoord.updateRhs(info, tmpView);
   Profiler::RegionStop<QUICC_DETAIL_PROF_LVL>("Timestep-input:updateRhs");

                  // If required set inhomogenous boundary condition value
                  if(cinfo.hasBoundaryValue())
                  {
                     // Set boundary value
                     std::visit(
                           [&](auto&& p)
                           {
                           throw std::logic_error("NOT YET IMPLEMENTED");
                           // Equations::setBoundaryValue(spEq->res(), cinfo, spEq->boundaryKernel(id.second), p->dom(0).perturbation(), id.second, (*solveIt)->rInhomogeneous(i), i, (*solveIt)->startRow(id,i));
                           }, eqIt->spUnknown());
                     //               this->mSolverCoord.updateInhomogeneous(info);
                  }
               }
            }
         }
      }
   };

   processInput(scalEq);
   processInput(vectEq);
}

template <typename TScheme>
void InterfaceViews<TScheme>::transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
{
   Profiler::RegionFixture<2> fix("Timestep-output");

   auto processOutput = [this](auto&& eq_range)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      for(auto& eqIt: make_range(eq_range))
      {
         for(auto& compIt: make_range(eqIt->spectralRange()))
         {
            // Get field identity
            myId = std::make_pair(eqIt->name(), compIt);

            const auto& cinfo = eqIt->couplingInfo(myId.second);

            if(eqIt->solveTiming() == SolveTiming::Prognostic::id())
            {
               DebuggerMacro_msg("Get timestepper solution for " + PhysicalNames::Coordinator::tag(myId.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(myId.second)) + ")", 6);

   Profiler::RegionStart<QUICC_DETAIL_PROF_LVL>("Timestep-output:setZero");
               // return zero
               for(std::size_t i = 0; i < static_cast<std::size_t>(cinfo.fieldStart()); i++)
               {
                  DecoupledZMatrix tmp(cinfo.galerkinN(i), cinfo.rhsCols(i));
                  tmp.setZero();

                  const SparseMatrix* pOp;
                  if(cinfo.isGalerkin())
                  {
                     pOp = &eqIt->galerkinStencil(myId.second, i);
                  }
                  std::visit(
                        [&](auto&& p)
                        {
                           Equations::storeSolution(p->rDom(0).rPerturbation(), eqIt->res(), cinfo, pOp, eqIt->solutionUpdater(myId.second), myId.second, tmp, i, 0);
                        }, eqIt->spUnknown());
               }
   Profiler::RegionStop<QUICC_DETAIL_PROF_LVL>("Timestep-output:setZero");

   Profiler::RegionStart<QUICC_DETAIL_PROF_LVL>("Timestep-output:allocateTemp");
               // Allocate temporary storage
               std::uint32_t mem_size = 0;
               for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
               {
                  mem_size = std::max(mem_size, static_cast<std::uint32_t>(cinfo.galerkinN(i))*static_cast<std::uint32_t>(cinfo.rhsCols(i)));
               }
               Memory::MemBlock<MHDComplex> data(mem_size, this->_mem.get());
   Profiler::RegionStop<QUICC_DETAIL_PROF_LVL>("Timestep-output:allocateTemp");

               for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
               {
   Profiler::RegionStart<QUICC_DETAIL_PROF_LVL>("Timestep-output:createInfo");
                  auto info = createInfo(cinfo, i);
   Profiler::RegionStop<QUICC_DETAIL_PROF_LVL>("Timestep-output:createInfo");

   Profiler::RegionStart<QUICC_DETAIL_PROF_LVL>("Timestep-output:getSolution");
                  std::uint32_t mem_rows = static_cast<std::uint32_t>(cinfo.galerkinN(i));
                  std::uint32_t mem_cols = static_cast<std::uint32_t>(cinfo.rhsCols(i));
                  using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
                  std::array<std::uint32_t, 2> dimensions {mem_rows, mem_cols};
                  View::View<MHDComplex, View::Attributes<dense2D>> tmpView(data, dimensions);
                  this->mSolverCoord.getSolution(tmpView, info);
   Profiler::RegionStop<QUICC_DETAIL_PROF_LVL>("Timestep-output:getSolution");

   Profiler::RegionStart<QUICC_DETAIL_PROF_LVL>("Timestep-output:copyView");
                  DecoupledZMatrix tmp(cinfo.galerkinN(i), cinfo.rhsCols(i));
                  Views::details::computeSet(tmp, tmpView, 0);
   Profiler::RegionStop<QUICC_DETAIL_PROF_LVL>("Timestep-output:copyView");

   Profiler::RegionStart<QUICC_DETAIL_PROF_LVL>("Timestep-output:storeSolution");
                  const SparseMatrix* pOp;
                  if(cinfo.isGalerkin())
                  {
                     pOp = &eqIt->galerkinStencil(myId.second, i);
                  }
                  std::visit(
                        [&](auto&& p)
                        {
                        Equations::storeSolution(p->rDom(0).rPerturbation(), eqIt->res(), cinfo, pOp, eqIt->solutionUpdater(myId.second), myId.second, tmp, i, 0);
                        }, eqIt->spUnknown());
   Profiler::RegionStop<QUICC_DETAIL_PROF_LVL>("Timestep-output:storeSolution");
               }

               // Apply constraint on solution
               auto changedSolution = eqIt->applyConstraint(myId.second, SolveTiming::After::id());

               // Update timestepper solver solution if constraint modified it
               if(changedSolution)
               {
                  auto corr_ = eqIt->correctionConstraint(myId.second, SolveTiming::After::id());
                  for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
                  {
                     auto info = createInfo(cinfo, i);

                     // Get effective corrections
                     auto corr = Equations::correctSolution(eqIt->res(), cinfo, myId.second, corr_, i, 0);

                     this->mSolverCoord.updateSolution(info, corr);
                  }
               }
            }
         }
      }
   };

   processOutput(scalEq);
   processOutput(vectEq);
}

template <typename TScheme>
void InterfaceViews<TScheme>::adaptTimestep(const Matrix& cfl)
{
   // Process CFL information
   this->processCfl(cfl, this->mSolverCoord.error(), this->mspScheme->order());

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
   this->writeCfl();
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
   // Timestep scheme
   std::stringstream oss;
   oss << "Timestepper: " << this->mspScheme->name() << " ("
       << this->mspScheme->order() << ")";

   this->printInfo(stream, oss.str());
}

// \todo  move to separate file
inline void buildTimestepMatrixWrapper(std::map<std::size_t, DecoupledZSparse>& ops, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp,
   const int idx)
{
   auto buildOp = [&](const std::size_t opId, const std::size_t bcId)
   {
      auto ret = ops.insert(std::make_pair(opId, DecoupledZSparse()));
      spEq->buildModelMatrix(ret.first->second, opId, comp, idx, bcId);
   };

   bool isSplit = spEq->couplingInfo(comp).isSplitEquation();

   using namespace ModelOperator;
   using namespace ModelOperatorBoundary;
   // Compute model's linear operator (without Tau lines)
   buildOp(ImplicitLinear::id(), SolverNoTau::id());

   // Compute model's time operator (without Tau lines)
   buildOp(Time::id(), SolverNoTau::id());

   // Compute model's quasi-inverse operator (without Tau lines)
   buildOp(QuasiInverse::id(), SolverNoTau::id());

   // Compute model's tau line boundary operator
   buildOp(Boundary::id(), SolverHasBc::id());

   // If equation was split into two lower order systems
   if (isSplit)
   {
      // Compute model's split linear operator (without Tau lines)
      buildOp(SplitImplicitLinear::id(), SolverNoTau::id());

      // Compute model's split quasi-inverse operator (without Tau lines)
      buildOp(SplitQuasiInverse::id(), SolverNoTau::id());

      // Compute model's tau line boundary operator for split operator
      buildOp(SplitBoundary::id(), SolverHasBc::id());

      // Compute model's tau line boundary value for split operator
      buildOp(SplitBoundaryValue::id(), SolverNoTau::id());
   }
}

} // namespace PredictorCorrector
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_PREDICTORCORRECTOR_INTERFACEVIEWS_HPP
