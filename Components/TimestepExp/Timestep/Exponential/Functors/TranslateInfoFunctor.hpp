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
#include "QuICC/Equations/CouplingInformation.hpp"
#include "Timestep/Exponential/TimestepperInfo.hpp"
#include "Timestep/Exponential/CreateInfo.hpp"
#include "Timestep/Exponential/Functors/BaseFunctor.hpp"
#include "Timestep/Exponential/Functors/FunctorData.hpp"
#include "Timestep/Exponential/BuildTimestepMatrix.hpp"
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
 * @brief Translate into timestepper info
 */
class TranslateInfoFunctor: public BaseFunctor
{
   public:
      TranslateInfoFunctor(std::vector<TimestepperInfo>& infos, std::shared_ptr<FunctorData> spData, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem) : BaseFunctor(idMap, mem), infos(infos), mspData(spData) {};
      ~TranslateInfoFunctor() = default;
      void operator()(const SpectralFieldId& id);
   private:
      std::vector<TimestepperInfo>& infos;
      std::shared_ptr<FunctorData> mspData;
};

inline void TranslateInfoFunctor::operator()(const SpectralFieldId& myId)
{
   assert(this->mspData->spRes);
   auto&& res = *this->mspData->spRes;
   auto&& bcIdMap = this->mspData->bcIdMap;
   auto&& eqParamsMap = this->mspData->eqParamsMap;
   assert(this->mspData->spBackend);
   auto&& backend = *this->mspData->spBackend;
   const auto& cinfo = this->mspData->cInfos.at(myId);

   DebuggerMacro_msg("Creating timesteppers for " + PhysicalNames::Coordinator::tag(myId.first) + "(" + Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(myId.second)) + ")", 2);

   auto& fId = *this->mpIdMap;
   fId.emplace(myId, fId.size());

   auto info = createInfo(cinfo, cinfo.fieldStart(), fId.at(myId));
   info.rows = 0;
   info.blockN = 0;

   for(std::size_t i = cinfo.fieldStart(); i < static_cast<std::size_t>(cinfo.nSystems()); i++)
   {
      auto t = createInfo(cinfo, i, fId.at(myId));

      // Set operators
      std::map<std::size_t, DecoupledZSparse> ops;
      buildTimestepMatrix(ops, myId.second, i, res, backend, cinfo, bcIdMap, eqParamsMap);

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
