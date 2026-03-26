/**
 * @file TranslateInfoFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <memory>

// Project includes
//
#include "Memory/MemoryResource.hpp"
#include "Timestep/Exponential/TimestepperInfo.hpp"
#include "Timestep/Exponential/CreateInfo.hpp"
#include "Timestep/Exponential/Functors/BaseFunctor.hpp"
#include "Timestep/Exponential/BuildTimestepMatrixWrapper.hpp"
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Translate into timestepper info
 */
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

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
