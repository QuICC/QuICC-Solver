#include "TestSuite/ViewMeta.hpp"
#include "TestSuite/Io.hpp"

namespace QuICC {
namespace TestSuite {

dimsAndMeta readDimsAndMeta(const std::string path, const std::string dist,
   const std::string id)
{
   int rank, ranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);

   dimsAndMeta ret;

   std::string base = path + dist + "/P_id" + id + "_np" +
                      std::to_string(ranks) + "_r" + std::to_string(rank);
   std::string fileName = base + "_stage0_meta.dat";
   std::vector<MHDFloat> dbJW;
   QuICC::TestSuite::readList(dbJW, fileName);
   std::uint32_t modsN = dbJW[0];
   std::uint32_t N = dbJW[1];

   fileName = base + "_stage1_meta.dat";
   std::vector<MHDFloat> dbAL;
   QuICC::TestSuite::readList(dbAL, fileName);
   std::uint32_t modsL = dbAL[0];
   std::uint32_t L = dbAL[1];

   fileName = base + "_stage2_meta.dat";
   std::vector<MHDFloat> dbFr;
   QuICC::TestSuite::readList(dbFr, fileName);
   std::uint32_t modsM = dbFr[0];
   std::uint32_t M = dbFr[1];

   // Grid dimensions
   constexpr std::uint32_t dim = 3;
   // RThetaPhi - NLM - v012
   ret.physDims = std::array<std::uint32_t, dim>{N, L, M};
   ret.modsDims = std::array<std::uint32_t, dim>{modsN, modsL, modsM};
   // Unpack metadata Fourier stage
   ret.metaFT = std::move(unPackMeta(dbFr, N));
   // Unpack metadata AL stage
   ret.metaAL = std::move(unPackMeta(dbAL, modsM));
   // Unpack metadata JW stage
   ret.metaJW = std::move(unPackMetaJW(dbJW, modsN));

   return ret;
}

View::ptrAndIdx unPackMeta(const std::vector<MHDFloat>& db,
   std::uint32_t maxLayers)
{
   // Unpack metadata
   std::uint32_t nLayers = db[2];
   // Pointers
   // note that index order is reversed
   std::vector<std::uint32_t> ptr(maxLayers + 1);
   ptr[0] = 0;
   for (std::size_t i = 1; i < ptr.size(); ++i)
   {
      ptr[i] += ptr[i - 1];
      for (std::uint32_t l = 0; l < nLayers; ++l)
      {
         if (db[3 + l * 2] == i - 1)
         {
            auto nCols = db[4 + l * 2];
            ptr[i] += nCols;
            continue;
         }
      }
   }
   // Indices
   auto totCols = ptr[maxLayers];
   std::vector<std::uint32_t> idx(totCols);
   std::uint32_t currCols = 0;
   for (std::size_t i = 0; i < ptr.size() - 1; ++i)
   {
      auto nCols = ptr[i + 1] - ptr[i];
      for (std::size_t c = 0; c < nCols; ++c)
      {
         idx[c + ptr[i]] = db[5 + (nLayers - 1) * 2 + (currCols) * 2 + c * 2];
      }
      currCols += nCols;
   }
   return {ptr, idx};
}


View::ptrAndIdx unPackMetaJW(const std::vector<MHDFloat>& db,
   std::uint32_t maxLayers)
{
   // Unpack metadata
   std::uint32_t nLayers = db[2];
   // Pointers
   // note that index order is reversed
   std::vector<std::uint32_t> ptr(maxLayers + 1);
   ptr[0] = 0;
   for (std::size_t i = 1; i < ptr.size(); ++i)
   {
      ptr[i] += ptr[i - 1];
      for (std::uint32_t l = 1; l <= nLayers; ++l)
      {
         if (db[3 + (nLayers - l) * 2] == maxLayers - i)
         {
            auto nCols = db[4 + (nLayers - l) * 2];
            ptr[i] += nCols;
            continue;
         }
      }
   }
   // Indices
   auto totCols = ptr[maxLayers];
   std::vector<std::uint32_t> idx(totCols);
   for (std::size_t i = 0; i < ptr.size() - 1; ++i)
   {
      auto nCols = ptr[i + 1] - ptr[i];
      for (std::size_t c = 0; c < nCols; ++c)
      {
         idx[c + ptr[i]] =
            db[5 + (nLayers - 1) * 2 + (totCols - 1) * 2 - (nCols - 1 - c) * 2];
      }
      totCols -= nCols;
   }
   return {ptr, idx};
}

} // namespace TestSuite
} // namespace QuICC
