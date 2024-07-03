#pragma once

#include <cstdint>

#include "Graph/Shims/MlirShims.hpp"
#include "Types/BasicTypes.hpp"


namespace QuICC
{
namespace Graph
{

ptrAndIdx unPackMeta(const std::vector<MHDFloat>& db, std::uint32_t maxLayers);

ptrAndIdx unPackMetaJW(const std::vector<MHDFloat>& db, std::uint32_t maxLayers);


} // namespace Graph
} // namespace QuICC

