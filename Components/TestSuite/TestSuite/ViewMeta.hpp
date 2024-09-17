/**
 * @file ViewMeta.hpp
 * @brief Utilities to setup views meta data for tests
 */
#pragma once

// System includes
//
#include <cstdint>
#include <vector>
#ifdef QUICC_MPI
#include <mpi.h>
#endif

// Project includes
//
#include "Types/BasicTypes.hpp"
#include "View/ViewUtils.hpp"


namespace QuICC {
namespace TestSuite {

constexpr std::uint32_t dim = 3;
/// @brief struct to pack promblem dimensions and metadata
struct dimsAndMeta
{
   /// @brief Physical dimensions
   /// RThetaPhi - v012
   std::array<std::uint32_t, dim> physDims;
   /// @brief Spectral dimensions
   /// NLM - v012
   std::array<std::uint32_t, dim> modsDims;
   View::ptrAndIdx metaFT;
   View::ptrAndIdx metaAL;
   View::ptrAndIdx metaJW;
   View::ptrAndIdx metaIM;
};

/// @brief util to setup problem dimensions and meta from files
/// @param path path to meta data files
/// @param dist parallelization scheme
/// @return
dimsAndMeta readDimsAndMeta(const std::string path, const std::string dist,
   const std::string id);

/// @brief util to reconstruc meta from file
/// @param db
/// @param maxLayers
/// @return
View::ptrAndIdx unPackMeta(const std::vector<MHDFloat>& db,
   std::uint32_t maxLayers);

} // namespace TestSuite
} // namespace QuICC
