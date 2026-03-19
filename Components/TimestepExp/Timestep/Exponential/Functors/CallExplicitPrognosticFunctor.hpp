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
#include "Timestep/Exponential/InterfaceFunctors.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

template <typename TTsFunc>
class CallExplicitPrognosticFunctor: public Pseudospectral::ExplicitPrognosticFunctor
{
   public:
      typedef std::map<SpectralFieldId, std::size_t> IdMap;
      CallExplicitPrognosticFunctor(std::shared_ptr<TTsFunc> tsFunc, const std::size_t regId, const std::size_t regCol, const int fixedIt, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem): tsFunc(tsFunc), regId(regId), regCol(regCol), fixedIt(fixedIt), mpIdMap(idMap), _mem(mem){};
      virtual ~CallExplicitPrognosticFunctor() = default;
      void operator()(const std::size_t opId, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, ScalarVariable_map& scalVar, VectorVariable_map& vectVar) final;
   protected:
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
   auto pvFunc = std::make_shared<OpFunctor>(this->tsFunc, opId, regId, regCol, scalVar, vectVar);
   auto pFunc = std::make_shared<LinearInputFunctor<OpFunctor>>(pvFunc, opId, this->mpIdMap, this->_mem);
   auto aFunc = std::make_shared<DoNothingFunctor>();
   ProcessRangeFunctor processor(bFunc, pFunc, aFunc, fixedIt);
   processor(scalEq);
   processor(vectEq);
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
