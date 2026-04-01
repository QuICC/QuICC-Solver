/**
 * @file InputFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <memory>

// Project includes
//
#include "Types/BasicTypes.hpp"
#include "View/Attributes.hpp"
#include "View/ViewDense.hpp"
#include "Memory/Memory.hpp"
#include "Memory/MemoryResource.hpp"
#include "Timestep/Exponential/CreateInfo.hpp"
#include "Timestep/Exponential/Functors/BaseFunctor.hpp"
#include "Timestep/Exponential/Functors/FunctorData.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#endif //QUICC_DEBUG

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Input functor
 */
template <typename TFunc>
class InputFunctor: public BaseFunctor
{
   public:
      InputFunctor(std::shared_ptr<FunctorData> spData, std::shared_ptr<TFunc> vFunc, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem) : BaseFunctor(idMap, mem), spData(spData), vFunc(vFunc){};
      virtual ~InputFunctor() = default;
      void operator()(const SpectralFieldId& id);
   protected:
      std::shared_ptr<FunctorData> spData; 
      std::shared_ptr<TFunc> vFunc;
};

template <typename TFunc>
void InputFunctor<TFunc>::operator()(const SpectralFieldId& myId)
{
   const auto& eqData = *spData;
   assert(eqData.cInfos.count(myId) > 0);
   const auto& cinfo = eqData.cInfos.at(myId);
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

      (*vFunc)(tmpView, myId, cinfo, info, i);
   }
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
