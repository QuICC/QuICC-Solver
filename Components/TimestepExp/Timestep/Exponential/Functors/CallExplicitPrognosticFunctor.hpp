/**
 * @file CallExplicitPrognosticFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <memory>

// Project includes
//
#include "Memory/MemoryResource.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Pseudospectral/Coordinator.hpp"
#include "Timestep/Exponential/Functors/DoNothingFunctor.hpp"
#include "Timestep/Exponential/Functors/FunctorData.hpp"
#include "Timestep/Exponential/Functors/GetLinearInputFunctor.hpp"
#include "Timestep/Exponential/Functors/LinearInputFunctor.hpp"
#include "Timestep/Exponential/Functors/ProcessRangeFunctor.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

template <typename TTsFunc>
class CallExplicitPrognosticFunctor: public Pseudospectral::ExplicitPrognosticFunctor
{
   public:
      typedef std::map<SpectralFieldId, std::size_t> IdMap;
      /**
       * @brief ctor
       */
      CallExplicitPrognosticFunctor(std::shared_ptr<FunctorData> spData, std::shared_ptr<TTsFunc> tsFunc, const std::size_t regId, const std::size_t regCol, const int fixedIt, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem): spData(spData), tsFunc(tsFunc), regId(regId), regCol(regCol), fixedIt(fixedIt), mpIdMap(idMap), _mem(mem){};

      /**
       * @brief ctor
       */
      virtual ~CallExplicitPrognosticFunctor() = default;

      void operator()(const std::size_t opId, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, ScalarVariable_map& scalVar, VectorVariable_map& vectVar) final;

   protected:
      std::shared_ptr<FunctorData> spData;
      std::shared_ptr<TTsFunc> tsFunc;
      const std::size_t regId;
      const std::size_t regCol;
      const int fixedIt;
   /**
    * @brief Shared field ID to solver field id
    */
   std::shared_ptr<IdMap> mpIdMap;

   /**
    * @brief
    */
   std::shared_ptr<Memory::memory_resource> _mem;
};

template <typename TTsFunc>
void CallExplicitPrognosticFunctor<TTsFunc>::operator()(const std::size_t opId, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, ScalarVariable_map& scalVar, VectorVariable_map& vectVar)
{
   Profiler::RegionFixture<2> fix("Timestep-explicitInput");

   auto bFunc = std::make_shared<DoNothingFunctor>();
   using OpFunctor = GetLinearInputFunctor<TTsFunc>;
   auto pvFunc = std::make_shared<OpFunctor>(this->spData, this->tsFunc, opId, regId, regCol, scalVar, vectVar);
   auto pFunc = std::make_shared<LinearInputFunctor<OpFunctor>>(this->spData, pvFunc, opId, this->mpIdMap, this->_mem);
   auto aFunc = std::make_shared<DoNothingFunctor>();
   ProcessRangeFunctor processor(bFunc, pFunc, aFunc, fixedIt);
   processor(this->spData->eqInfos);
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
