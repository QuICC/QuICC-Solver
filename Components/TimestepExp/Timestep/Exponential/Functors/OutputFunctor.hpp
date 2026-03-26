/**
 * @file OutputFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "Types/Typedefs.hpp"
#include "View/Attributes.hpp"
#include "View/ViewDense.hpp"
#include "Memory/Memory.hpp"
#include "Memory/MemoryResource.hpp"
#include "Timestep/Exponential/Functors/BaseFunctor.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Output functor
 */
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

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
