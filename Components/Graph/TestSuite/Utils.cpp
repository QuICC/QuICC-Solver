#include "Utils.hpp"

namespace QuICC
{
namespace Graph
{

ptrAndIdx unPackMeta(const std::vector<MHDFloat>& db, std::uint32_t maxLayers)
{
  // Unpack metadata
  std::uint32_t nLayers = db[2];
  // Pointers
  // note that index order is reversed
  std::vector<std::uint32_t> ptr(maxLayers+1);
  ptr[0] = 0;
  for (std::size_t i = 1; i < ptr.size(); ++i)
  {
    ptr[i] += ptr[i-1];
    for (std::uint32_t l = 0; l < nLayers; ++l)
    {
      if (db[3 + l*2] == i - 1)
      {
        auto nCols = db[4 + l*2];
        ptr[i] += nCols;
        continue;
      }
    }
  }
  // Indices
  auto totCols = ptr[maxLayers];
  std::vector<std::uint32_t> idx(totCols);
  std::uint32_t currCols = 0;
  for (std::size_t i = 0; i < ptr.size()-1; ++i)
  {
    auto nCols = ptr[i+1] - ptr[i];
    for (std::size_t c = 0; c < nCols; ++c)
    {
      idx[c+ptr[i]] = db[5 + (nLayers-1)*2 + (currCols)*2 +c*2];
    }
    currCols += nCols;
  }
  return {ptr, idx};
}


ptrAndIdx unPackMetaJW(const std::vector<MHDFloat>& db, std::uint32_t maxLayers)
{
  // Unpack metadata
  std::uint32_t nLayers = db[2];
  // Pointers
  // note that index order is reversed
  std::vector<std::uint32_t> ptr(maxLayers+1);
  ptr[0] = 0;
  for (std::size_t i = 1; i < ptr.size(); ++i)
  {
    ptr[i] += ptr[i-1];
    for (std::uint32_t l = 1; l <= nLayers; ++l)
    {
      if (db[3 + (nLayers-l)*2] == maxLayers - i)
      {
        auto nCols = db[4 + (nLayers-l)*2];
        ptr[i] += nCols;
        continue;
      }
    }
  }
  // Indices
  auto totCols = ptr[maxLayers];
  std::vector<std::uint32_t> idx(totCols);
  for (std::size_t i = 0; i < ptr.size()-1; ++i)
  {
    auto nCols = ptr[i+1] - ptr[i];
    for (std::size_t c = 0; c < nCols; ++c)
    {
      idx[c+ptr[i]] = db[5 + (nLayers-1)*2 + (totCols-1)*2 - (nCols-1-c)*2];
    }
    totCols -= nCols;
  }
  return {ptr, idx};
}

} // namespace Graph
} // namespace QuICC
