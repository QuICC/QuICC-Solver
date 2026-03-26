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
      LinearInputFunctor(std::shared_ptr<TFunc> vFunc, const std::size_t opId, std::shared_ptr<BaseFunctor::IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem) : InputFunctor<TFunc>(vFunc, idMap, mem), opId(opId){};
      virtual ~LinearInputFunctor() = default;
      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt);
   protected:
      const std::size_t opId;
};

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


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
