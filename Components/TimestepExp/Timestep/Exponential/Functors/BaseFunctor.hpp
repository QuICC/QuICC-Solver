/**
 * @file BaseFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <memory>
#include <map>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "Memory/MemoryResource.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Base timestepping functors
 */
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

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
