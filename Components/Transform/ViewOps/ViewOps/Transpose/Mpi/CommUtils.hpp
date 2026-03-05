/**
 * @file Comm.hpp
 * @brief Methods for mpi enabled transform
 */
#pragma once

// External includes
//
#include <array>
#include <cassert>
#include <memory>
#include <mpi.h>
#include <vector>

// Project includes
//
#include "Environment/MpiTypes.hpp"
#include "Memory/Cpu/NewDelete.hpp"
#include "Memory/Memory.hpp"
#include "View/View.hpp"
#include "ViewOps/Transpose/Mpi/Tags.hpp"
#include "ViewOps/Transpose/Packing.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "Cuda/CudaUtil.hpp"
#include "Memory/Cuda/Malloc.hpp"
#endif

namespace QuICC {
namespace Transpose {
namespace Mpi {

/// \todo move Mpi utils out of Enviroment and unify
namespace details {
inline void mpiAssert(int ierr)
{
#ifndef NDEBUG
   if (ierr != MPI_SUCCESS)
   {
      throw std::runtime_error("Mpi failed.");
   }
#endif
}
} // namespace details


/// @brief point coordinate dimensions
constexpr int dimSize = 3;

/// @brief point coordinate type
using point_t = std::array<int, dimSize>;

/// @brief Build send or recv displacement
/// @param absCooNew ending coordinates
/// @param absCooOld starting coordinates
/// @param comm mpi communicator spanning all involved ranks
/// @return send/recv displacement
std::vector<std::vector<int>> getDispls(const std::vector<point_t>& absCooNew,
   const std::vector<point_t>& absCooOld, const MPI_Comm comm = MPI_COMM_WORLD);


/// @brief Build set of communicating ranks from displacements
/// This is part of the alltoallw setup
/// @param sendDispls sending displacements
/// @param recvDispls receiving displacements
/// @return set of communicating ranks
std::vector<int> getReducedRanksSet(
   const std::vector<std::vector<int>>& sendDispls,
   const std::vector<std::vector<int>>& recvDispls,
   const MPI_Comm comm = MPI_COMM_WORLD);

/// @brief Reduce diplacements to the to the reduced set
/// @param sendDispls sending displacements
/// @param recvDispls receiving displacements
/// @param redSet reduced set of ranks
void redDisplsFromSet(std::vector<std::vector<int>>& sendDispls,
   std::vector<std::vector<int>>& recvDispls, const std::vector<int>& redSet);

/// @brief Get sub communicator from the reduced set
/// @param redSet reduced set of ranks
/// @param comm all world communicator
/// @return sub communicator
MPI_Comm getSubComm(const std::vector<int>& redSet,
   const MPI_Comm comm = MPI_COMM_WORLD);

/// @brief Get count vector based on displacements
/// @param displs displacements
/// @return count vector
template <class TAG>
std::vector<int> getCount(const std::vector<std::vector<int>>& displs)
{
   std::vector<int> count(displs.size());
   for (std::size_t i = 0; i < displs.size(); ++i)
   {
      count[i] = displs[i].size();
   }
   return count;
}

} // namespace Mpi
} // namespace Transpose
} // namespace QuICC
