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
#include "QuICC/Diagnostics/Coordinator.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "Timestep/Exponential/TimestepperInfo.hpp"
#include "Timestep/Exponential/CreateInfo.hpp"
#include "Timestep/Exponential/Functors/BaseFunctor.hpp"
#include "Timestep/Exponential/BuildTimestepMatrix.hpp"
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
      typedef std::map<SpectralFieldId, Equations::CouplingInformation> CInfoMap;
      TranslateInfoFunctor(std::vector<TimestepperInfo>& infos, std::shared_ptr<CInfoMap> cinfoMap, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem) : BaseFunctor(idMap, mem), mpCInfoMap(cinfoMap), infos(infos), mspRes(nullptr), mspBcIdMap(nullptr), mspEqParamsMap(nullptr), mspBackend(nullptr) {};
      ~TranslateInfoFunctor() = default;
      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt);
      std::shared_ptr<Resolution> spRes() const;
      std::shared_ptr<std::map<std::size_t, std::size_t>> spBcIdMap() const;
      std::shared_ptr<std::map<std::size_t, NonDimensional::SharedINumber>> spEqParamsMap() const;
      std::shared_ptr<Model::IModelBackend> spBackend() const;
   private:
      std::shared_ptr<CInfoMap> mpCInfoMap;
      std::vector<TimestepperInfo>& infos;
      std::shared_ptr<Resolution> mspRes;
      std::shared_ptr<std::map<std::size_t, std::size_t>> mspBcIdMap;
      std::shared_ptr<std::map<std::size_t, NonDimensional::SharedINumber>> mspEqParamsMap;
      std::shared_ptr<Model::IModelBackend> mspBackend;
};

template <typename TEqIt>
void TranslateInfoFunctor::operator()(const SpectralFieldId& myId, TEqIt& eqIt)
{
   auto& fId = *this->mpIdMap;
   auto& ciMap = *this->mpCInfoMap;

   // Set Resolution
   if(this->mspRes== nullptr)
   {
      this->mspRes = eqIt->spRes();
   }
   auto&& res = *this->mspRes;

   // Set BcIdMap
   if(this->mspBcIdMap== nullptr)
   {
      this->mspBcIdMap = std::make_shared<std::map<std::size_t, std::size_t>>();
      *this->mspBcIdMap = eqIt->bcIds().map();
   }
   auto&& bcIdMap = *this->mspBcIdMap;

   // Set EqParamsMap
   if(this->mspEqParamsMap== nullptr)
   {
      this->mspEqParamsMap = std::make_shared<std::map<std::size_t, NonDimensional::SharedINumber>>();
      *this->mspEqParamsMap = eqIt->eqParams().map();
   }
   auto&& eqParamsMap = *this->mspEqParamsMap;

   // Set backend
   if(this->mspBackend == nullptr)
   {
      this->mspBackend = eqIt->spBackend();
   }
   auto&& backend = *this->mspBackend;

   ciMap.emplace(myId, eqIt->couplingInfo(myId.second));
   const auto& cinfo = ciMap.at(myId);
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

inline std::shared_ptr<Model::IModelBackend> TranslateInfoFunctor::spBackend() const
{
   return this->mspBackend;
}

inline std::shared_ptr<Resolution> TranslateInfoFunctor::spRes() const
{
   return this->mspRes;
}

inline std::shared_ptr<std::map<std::size_t, std::size_t>> TranslateInfoFunctor::spBcIdMap() const
{
   return this->mspBcIdMap;
}

inline std::shared_ptr<std::map<std::size_t, NonDimensional::SharedINumber>> TranslateInfoFunctor::spEqParamsMap() const
{
   return this->mspEqParamsMap;
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
