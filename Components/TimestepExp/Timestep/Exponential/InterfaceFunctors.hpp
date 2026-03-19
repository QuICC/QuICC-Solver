/**
 * @file InterfaceFunctors.hpp
 * @brief Implementation of functors used in Interface
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_INTERFACEFUNCTORS_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_INTERFACEFUNCTORS_HPP

// System includes
//
#include <memory>
#include <map>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Memory/MemoryResource.hpp"
#include "View/Attributes.hpp"
#include "View/ViewDense.hpp"
#include "Memory/Memory.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "Timestep/Exponential/TimestepperInfo.hpp"
#include "QuICC/SolveTiming/Before.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "Timestep/Exponential/details/TimesteppperTools.hpp"
#include "QuICC/Equations/CopyUnknown.hpp"
#include "QuICC/Equations/CorrectSolution.hpp"
#include "QuICC/Equations/AddSource.hpp"
#include "QuICC/Equations/SolveStencilUnknown.hpp"
#include "QuICC/Equations/ExplicitTerm.hpp"
#include "QuICC/IteratorRange.hpp"
#include "QuICC/Timestep/Interface.hpp"
#include "QuICC/ModelOperator/QuasiInverse.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoBc.hpp"
#include "QuICC/Tag/Operator/Qi.hpp"
#include "QuICC/Tag/Operator/Lhs.hpp"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#endif

namespace QuICC {

namespace Timestep {

namespace Exponential {

TimestepperInfo createInfo(const Equations::CouplingInformation& cinfo, const std::size_t idx, const std::size_t fieldIndex);

/**
 * @brief Wrapper to build timestepping matrices
 */
void buildTimestepMatrixWrapper(std::map<std::size_t, DecoupledZSparse>& ops, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

class DoNothingFunctor
{
   public:
      DoNothingFunctor() = default;
      ~DoNothingFunctor() = default;
      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt){};
      
};

template <typename TCoord>
class GetStepperFunctor
{
   public:
      typedef std::pair<typename TCoord::TimestepperType*, std::size_t> ReturnType;

      GetStepperFunctor(TCoord& coord): coord(coord) {};
      ~GetStepperFunctor() = default;
      ReturnType operator()(const TimestepperInfo& info);
   private:
      TCoord& coord;
};

template <typename TStepper>
class StepperWrapperFunctor
{
   public:
      typedef std::pair<TStepper*, std::size_t> ReturnType;

      StepperWrapperFunctor(const std::vector<std::size_t>& startArr): startArr(startArr), pStepper(nullptr){};
      ~StepperWrapperFunctor() = default;
      void setStepper(TStepper& ts);
      ReturnType operator()(const TimestepperInfo& info);
   private:
      std::vector<std::size_t> startArr;
      TStepper* pStepper;
};

class ApplyConstraintFunctor
{
   public:
      ApplyConstraintFunctor(const std::size_t t): timing(t) {};
      ~ApplyConstraintFunctor() = default;
      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt);
   private:
      std::size_t timing;
};

class BaseFunctor
{
public:
   /// Typedef for Field ID to solver field ID
   typedef std::map<SpectralFieldId, std::size_t> IdMap;
   BaseFunctor(std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem): mpIdMap(idMap), _mem(mem){};
      virtual ~BaseFunctor() = default;

protected:
   /**
    * @brief Shared field ID to solver field id
    */
   std::shared_ptr<IdMap> mpIdMap;

   /**
    * @brief
    */
   std::shared_ptr<Memory::memory_resource> _mem;
};

class TranslateInfoFunctor: public BaseFunctor
{
   public:
      TranslateInfoFunctor(std::vector<TimestepperInfo>& infos, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem) : BaseFunctor(idMap, mem), infos(infos) {};
      ~TranslateInfoFunctor() = default;
      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt);
   private:
      std::vector<TimestepperInfo>& infos;
};

template <typename TFunc>
class InputFunctor: public BaseFunctor
{
   public:
      InputFunctor(std::shared_ptr<TFunc> vFunc, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem) : BaseFunctor(idMap, mem), vFunc(vFunc){};
      virtual ~InputFunctor() = default;
      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt);
   protected:
      std::shared_ptr<TFunc> vFunc;
};

template <typename TFunc>
class LinearInputFunctor: public InputFunctor<TFunc>
{
   public:
      LinearInputFunctor(std::shared_ptr<TFunc> vFunc, const std::size_t opId, std::shared_ptr<BaseFunctor::IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem) : InputFunctor<TFunc>(vFunc, idMap, mem), opId(opId){};
      virtual ~LinearInputFunctor() = default;
      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt);
   protected:
      const std::size_t opId;
};

template <typename TFunc, typename TCorrFunc>
class OutputFunctor: public BaseFunctor
{
   public:
      OutputFunctor(std::shared_ptr<TFunc> vFunc, std::shared_ptr<TCorrFunc> cFunc, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem) : BaseFunctor(idMap, mem), vFunc(vFunc), cFunc(cFunc){};
      virtual ~OutputFunctor() = default;
      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt);
   protected:
      std::shared_ptr<TFunc> vFunc;
      std::shared_ptr<TCorrFunc> cFunc;
};

template <typename TTsFunc>
class InitSolutionFunctor
{
   public:
      using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
      typedef View::View<MHDComplex, View::Attributes<dense2D>> ViewType;

      InitSolutionFunctor(std::shared_ptr<TTsFunc> tsFunc): tsFunc(tsFunc){};
      ~InitSolutionFunctor() = default;
      template <typename TEqIt>
      void operator()(ViewType tmpView, const SpectralFieldId& id, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i);
   protected:
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
class GetInputFunctor
{
   public:
      using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
      typedef View::View<MHDComplex, View::Attributes<dense2D>> ViewType;

      GetInputFunctor(std::shared_ptr<TTsFunc> tsFunc, const std::size_t regId, const std::size_t col): regId(regId), col(col), tsFunc(tsFunc){};
      ~GetInputFunctor() = default;
      template <typename TEqIt>
      void operator()(ViewType tmpView, const SpectralFieldId& id, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i);

   protected:
      const std::size_t regId;
      const std::size_t col;
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
class GetLinearInputFunctor
{
   public:
      using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
      typedef View::View<MHDComplex, View::Attributes<dense2D>> ViewType;

      GetLinearInputFunctor(std::shared_ptr<TTsFunc> tsFunc, const std::size_t opId, const std::size_t regId, const std::size_t col, const Timestep::Interface::ScalarVariable_map& scalVar, const Timestep::Interface::VectorVariable_map& vectVar) : opId(opId), regId(regId), col(col), scalVar(scalVar), vectVar(vectVar), tsFunc(tsFunc){};
      ~GetLinearInputFunctor() = default;
      template <typename TEqIt>
      void operator()(ViewType tmpView, const SpectralFieldId& id, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i);
   private:
      const std::size_t opId;
      const std::size_t regId;
      const std::size_t col;
      const Timestep::Interface::ScalarVariable_map& scalVar;
      const Timestep::Interface::VectorVariable_map& vectVar;
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
class TransferOutputFunctor
{
   public:
      using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
      typedef View::View<MHDComplex, View::Attributes<dense2D>> ViewType;

      TransferOutputFunctor(std::shared_ptr<TTsFunc> tsFunc, const std::size_t regId, const std::size_t col) : regId(regId), col(col), tsFunc(tsFunc){};
      ~TransferOutputFunctor() = default;
      template <typename TEqIt>
      void operator()(ViewType tmpView, const SpectralFieldId& id, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i);

   protected:
      const std::size_t regId;
      const std::size_t col;
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TTsFunc>
class TransferCorrectionFunctor
{
   public:
      TransferCorrectionFunctor(std::shared_ptr<TTsFunc> tsFunc, const std::size_t regId, const std::size_t col): regId(regId), col(col), tsFunc(tsFunc){};
      ~TransferCorrectionFunctor() = default;
      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const BaseFunctor::IdMap& idMap);
   protected:
      const std::size_t regId;
      const std::size_t col;
      std::shared_ptr<TTsFunc> tsFunc;
};

template <typename TBefore, typename TPrognostic, typename TAfter>
class ProcessRangeFunctor
{
   public:
      ProcessRangeFunctor(std::shared_ptr<TBefore> bFunc, std::shared_ptr<TPrognostic> pFunc, std::shared_ptr<TAfter> aFunc, const int fixedIt): _bFunc(bFunc), _pFunc(pFunc), _aFunc(aFunc), mRestricted(fixedIt != -1), mFixedIt(fixedIt) {};
      ~ProcessRangeFunctor() = default;
      template <typename TRange>
      void operator()(const TRange& eq_range);
   private:
      std::shared_ptr<TBefore> _bFunc;
      std::shared_ptr<TPrognostic> _pFunc;
      std::shared_ptr<TAfter> _aFunc;
      const bool mRestricted;
      const int mFixedIt;
};

template <typename TBefore, typename TPrognostic, typename TAfter>
template <typename TRange>
void ProcessRangeFunctor<TBefore,TPrognostic,TAfter>::operator()(const TRange& eq_range)
{
   TBefore& bFunc = *_bFunc;
   TPrognostic& pFunc = *_pFunc;
   TAfter& aFunc = *_aFunc;

   // Storage for information and identity
   SpectralFieldId myId;

   // Loop over equation range
   for(auto& eqIt: make_range(eq_range))
   {
      if(this->mRestricted && eqIt->options().it() != this->mFixedIt)
      {
         continue;
      }

      // Loop over spectral components
      for(auto& compId: make_range(eqIt->spectralRange()))
      {
         // Get field identity
         myId = std::make_pair(eqIt->name(), compId);

         // Process before prognostic equation
         bFunc(myId, eqIt);

         // Process prognostic equation
         if(eqIt->solveTiming() == SolveTiming::Prognostic::id())
         {
            pFunc(myId, eqIt);
         }

         // Process after prognostic equation
         aFunc(myId, eqIt);
      }
   }
}

template <typename TCoord>
typename GetStepperFunctor<TCoord>::ReturnType GetStepperFunctor<TCoord>::operator()(const TimestepperInfo& info)
{
   return coord.getStepper(info);
}

template <typename TStepper>
typename StepperWrapperFunctor<TStepper>::ReturnType StepperWrapperFunctor<TStepper>::operator()(const TimestepperInfo& info)
{
   std::size_t start = startArr.at(info.fieldIndex) + info.matStart;
   auto ts = std::make_pair(pStepper, start);
   return ts;
}

template <typename TStepper>
void StepperWrapperFunctor<TStepper>::setStepper(TStepper& ts)
{
   this->pStepper = &ts;
}

template <typename TEqIt>
void ApplyConstraintFunctor::operator()(const SpectralFieldId& myId, TEqIt& eqIt)
{
   // Apply constraint on solution
   eqIt->applyConstraint(myId.second, timing);
}

template <typename TEqIt>
void TranslateInfoFunctor::operator()(const SpectralFieldId& myId, TEqIt& eqIt)
{
   auto& fId = *this->mpIdMap;

   const auto& cinfo = eqIt->couplingInfo(myId.second);
   DebuggerMacro_msg("Creating timesteppers for " + PhysicalNames::Coordinator::tag(myId.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(myId.second)) + ")", 2);

   fId.emplace(myId, fId.size());

   auto info = createInfo(cinfo, cinfo.fieldStart(), fId.at(myId));
   info.rows = 0;
   info.blockN = 0;

   for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
   {
      auto t = createInfo(cinfo, i, fId.at(myId));

      // Set operators
      std::map<std::size_t, DecoupledZSparse> ops;
      buildTimestepMatrixWrapper(ops, eqIt, myId.second, i);

      for(auto&& op: ops)
      {
         if(info.ops.count(op.first) == 0)
         {
            info.ops.emplace(op.first, std::map<std::size_t, std::pair<int, DecoupledZSparse>>());
         }

         std::size_t matId = info.blockN;
         info.ops.at(op.first).emplace(matId, std::make_pair(cinfo.galerkinN(i), op.second));
      }

      info.rows += t.rows;
      info.matIds.push_back(info.blockN);
      info.blockN += t.blockN;
   }
   infos.push_back(info);
}

template <typename TFunc>
template <typename TEqIt>
void InputFunctor<TFunc>::operator()(const SpectralFieldId& myId, TEqIt& eqIt)
{
   const auto& cinfo = eqIt->couplingInfo(myId.second);
   const auto& idMap = *this->mpIdMap;

   // Allocate temporary storage
   std::uint32_t mem_size = 0;
   for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
   {
      mem_size = std::max(mem_size, static_cast<std::uint32_t>(cinfo.galerkinN(i))*static_cast<std::uint32_t>(cinfo.rhsCols(i)));
   }
   Memory::MemBlock<MHDComplex> data(mem_size, this->_mem.get());

   DebuggerMacro_msg("Get timestepper input for " + PhysicalNames::Coordinator::tag(myId.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(myId.second)) + ")", 6);

   // Get timestep input
   std::size_t matStart = 0;
   for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
   {
      auto info = createInfo(cinfo, i, idMap.at(myId));
      info.matStart = matStart;
      matStart += info.blockN;

      std::uint32_t mem_rows = static_cast<std::uint32_t>(cinfo.galerkinN(i));
      std::uint32_t mem_cols = static_cast<std::uint32_t>(cinfo.rhsCols(i));
      using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
      std::array<std::uint32_t, 2> dimensions {mem_rows, mem_cols};
      View::View<MHDComplex, View::Attributes<dense2D>> tmpView(data, dimensions);

      (*vFunc)(tmpView, myId, eqIt, cinfo, info, i);
   }
}

template <typename TFunc>
template <typename TEqIt>
void LinearInputFunctor<TFunc>::operator()(const SpectralFieldId& myId, TEqIt& eqIt)
{
   const auto& cinfo = eqIt->couplingInfo(myId.second);

   DebuggerMacro_msg("Get linear timestepper input for " + PhysicalNames::Coordinator::tag(myId.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(myId.second)) + ")", 6);

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
      InputFunctor<TFunc>::operator()(myId, eqIt);
   }
}

template <typename TFunc, typename TCorrFunc>
template <typename TEqIt>
void OutputFunctor<TFunc, TCorrFunc>::operator()(const SpectralFieldId& myId, TEqIt& eqIt)
{
   const auto& cinfo = eqIt->couplingInfo(myId.second);
   const auto& idMap = *this->mpIdMap;

   DebuggerMacro_msg("Get timestepper solution for " + PhysicalNames::Coordinator::tag(myId.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(myId.second)) + ")", 6);

   // return zero
   for(std::size_t i = 0; i < static_cast<std::size_t>(cinfo.fieldStart()); i++)
   {
      DecoupledZMatrix tmp(cinfo.galerkinN(i), cinfo.rhsCols(i));
      tmp.setZero();

      eqIt->storeSolution(myId.second, tmp, i, 0);
   }

   // Allocate temporary storage
   std::uint32_t mem_size = 0;
   for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
   {
      mem_size = std::max(mem_size, static_cast<std::uint32_t>(cinfo.galerkinN(i))*static_cast<std::uint32_t>(cinfo.rhsCols(i)));
   }
   Memory::MemBlock<MHDComplex> data(mem_size, this->_mem.get());

   std::size_t matStart = 0;
   for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
   {
      auto info = createInfo(cinfo, i, idMap.at(myId));
      info.matStart = matStart;
      matStart += info.blockN;

      std::uint32_t mem_rows = static_cast<std::uint32_t>(cinfo.galerkinN(i));
      std::uint32_t mem_cols = static_cast<std::uint32_t>(cinfo.rhsCols(i));
      using dense2D = View::DimLevelType<View::dense_t, View::dense_t>;
      std::array<std::uint32_t, 2> dimensions {mem_rows, mem_cols};
      View::View<MHDComplex, View::Attributes<dense2D>> tmpView(data, dimensions);

      (*vFunc)(tmpView, myId, eqIt, cinfo, info, i);
   }

   // Feedback for correcting timestepper solutions
   (*cFunc)(myId, eqIt, cinfo, idMap);
}

template <typename TTsFunc>
template <typename TEqIt>
void InitSolutionFunctor<TTsFunc>::operator()(ViewType tmpView, const SpectralFieldId& myId, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i)
{
   if(cinfo.isGalerkin())
   {
      Equations::solveStencilUnknown(*eqIt, myId.second, tmpView, i, 0);
   }
   else
   {
      std::visit(
            [&](auto&& p)
            {
            Equations::copyUnknown(*eqIt, p->dom(0).perturbation(), myId.second, tmpView, i, 0, true, true, true);
            }, eqIt->spUnknown());
   }

   auto tsData = (*tsFunc)(info);
   auto&& pStepper = tsData.first;
   auto&& start = tsData.second;
   pStepper->setSolution(tmpView, start);
   pStepper->updateSolutions();
}

template <typename TTsFunc>
template <typename TEqIt>
void GetInputFunctor<TTsFunc>::operator()(ViewType tmpView, const SpectralFieldId& myId, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i)
{
   // Copy field values into timestepper input
   if(cinfo.hasNonlinear())
   {
      std::visit(
            [&](auto&& p)
            {
            Equations::copyUnknown(*eqIt, p->dom(0).perturbation(), myId.second, tmpView, i, 0, true, true, false);
            }, eqIt->spUnknown());
   }

   // Add source term
   std::visit(
         [&](auto&& p)
         {
         Equations::addSource(*eqIt, p->dom(0).perturbation(), myId.second, tmpView, i, 0);
         }, eqIt->spUnknown());

   // Add value to RHS
   auto tsData = (*tsFunc)(info);
   auto&& pStepper = tsData.first;
   auto&& start = tsData.second;
   pStepper->addQiData(tmpView, start, this->regId, this->col);

   // Enforce BC
   const auto& rows = tmpView.dims()[0];
   const auto& cols = tmpView.dims()[1];
   pStepper->enforceBoundaryConditions(rows, cols, start, this->regId, this->col);

   // If required set inhomogenous boundary condition value
   if(cinfo.hasBoundaryValue())
   {
      throw std::logic_error("Inhomogeneous boundary conditions are not supported!");
   }
}

template <typename TTsFunc>
template <typename TEqIt>
void GetLinearInputFunctor<TTsFunc>::operator()(ViewType tmpView, const SpectralFieldId& myId, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i)
{
   // Copy field values into timestepper input
   DecoupledZMatrix tmp(cinfo.tauN(i), cinfo.rhsCols(i));
   tmp.setZero();

   // Build range of operator
   auto r = make_range(cinfo.explicitRange(opId));

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

   details::computeSet(tmpView, tmp, cinfo.tauN(i) - cinfo.galerkinN(i));

   auto tsData = (*tsFunc)(info);
   auto&& pStepper = tsData.first;
   auto&& start = tsData.second;
   pStepper->addData(tmpView, start, this->regId, this->col);
}

template <typename TTsFunc>
template <typename TEqIt>
void TransferOutputFunctor<TTsFunc>::operator()(ViewType tmpView, const SpectralFieldId& myId, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const TimestepperInfo& info, const std::size_t i)
{
   auto tsData = (*tsFunc)(info);
   auto&& pStepper = tsData.first;
   auto&& start = tsData.second;
   pStepper->getData(tmpView, start, this->regId, this->col);

   DecoupledZMatrix tmp(cinfo.galerkinN(i), cinfo.rhsCols(i));
   details::computeSet(tmp, tmpView, 0);

   eqIt->storeSolution(myId.second, tmp, i, 0);
}

template <typename TTsFunc>
template <typename TEqIt>
void TransferCorrectionFunctor<TTsFunc>::operator()(const SpectralFieldId& myId, TEqIt& eqIt, const Equations::CouplingInformation& cinfo, const BaseFunctor::IdMap& idMap)
{
   // Apply constraint on solution
   auto changedSolution = eqIt->applyConstraint(myId.second, SolveTiming::After::id());

   // Update timestepper solver solution if constraint modified it
   if(changedSolution)
   {
      auto corr_ = eqIt->correctionConstraint(myId.second, SolveTiming::After::id());
      std::size_t matStart = 0;
      for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
      {
         auto info = createInfo(cinfo, i, idMap.at(myId));
         info.matStart = matStart;
         matStart += info.blockN;

         // Get effective corrections
         auto corr = Equations::correctSolution(*eqIt, myId.second, corr_, i, 0);

         auto tsData = (*tsFunc)(info);
         auto&& pStepper = tsData.first;
         auto&& start = tsData.second;
         pStepper->correctData(corr, cinfo.galerkinN(i), cinfo.rhsCols(i), start, this->regId, this->col);
         pStepper->updateSolutions();
      }
   }
}

inline TimestepperInfo createInfo(const Equations::CouplingInformation& cinfo, const std::size_t idx, const std::size_t fieldIndex)
{
   TimestepperInfo info;
   info.isComplex = false;
   info.fieldIndex = fieldIndex;
   info.solverIndex = 0;
   info.rows = 2*cinfo.systemN(idx)*cinfo.rhsCols(idx);
   info.cols = 1;
   info.blockN = 2*cinfo.galerkinN(idx)*cinfo.rhsCols(idx);


   return info;
}

inline void buildTimestepMatrixWrapper(std::map<std::size_t, DecoupledZSparse>& ops, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp,
   const int idx)
{
   auto buildOp = [&](const std::size_t opId, const std::size_t tId, const std::size_t bcId)
   {
      auto ret = ops.insert(std::make_pair(tId, DecoupledZSparse()));
      spEq->buildModelMatrix(ret.first->second, opId, comp, idx, bcId);
   };

   using namespace ModelOperator;
   using namespace ModelOperatorBoundary;
   // Compute model's quasi-inverse operator (with boundary conditions)
   buildOp(Time::id(), Tag::Operator::Lhs::id(), SolverNoTau::id());

   // Compute model's quasi-inverse operator (without boundary conditions)
   buildOp(QuasiInverse::id(), Tag::Operator::Qi::id(), SolverNoBc::id());

}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_INTERFACEFUNCTORS_HPP
