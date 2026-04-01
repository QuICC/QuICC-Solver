/**
 * @file LinearInputFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <memory>

// Project includes
//
#include "Memory/MemoryResource.hpp"
#include "Timestep/Exponential/Functors/InputFunctor.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Linear input functor
 */
template <typename TFunc>
class LinearInputFunctor: public InputFunctor<TFunc>
{
   public:
      LinearInputFunctor(std::shared_ptr<FunctorData> spData, std::shared_ptr<TFunc> vFunc, const std::size_t opId, std::shared_ptr<BaseFunctor::IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem) : InputFunctor<TFunc>(spData, vFunc, idMap, mem), opId(opId){};
      virtual ~LinearInputFunctor() = default;
      void operator()(const SpectralFieldId& id);
   protected:
      const std::size_t opId;
};

template <typename TFunc>
void LinearInputFunctor<TFunc>::operator()(const SpectralFieldId& myId)
{
   const auto& data = *this->spData;
   assert(data.cInfos.count(myId) > 0);
   const auto& cinfo = data.cInfos.at(myId);

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
      InputFunctor<TFunc>::operator()(myId);
   }
}


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
