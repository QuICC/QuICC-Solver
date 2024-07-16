/**
 * @file MpiUtils.hpp
 * @brief Methods for mpi enabled transform
 */
#pragma once

// External includes
//
#include <array>
#include <mpi.h>
#include <vector>

// Project includes
//

namespace QuICC {
namespace Transpose {

/// @brief point coordinate dimensions
constexpr int dimSize = 3;

/// @brief point coordinate type
using point_t = std::array<int, dimSize>;

/// @brief Container for Mpi communicator and types
/// @tparam TAG
template<class TAG>
class Comm
{
public:
    Comm(MPI_Comm comm = MPI_COMM_WORLD):_comm(comm){};

    /// @brief release resources
    /// like Types and Displacements
    ~Comm();

    void setComm(std::vector<point_t>& cooNew, const std::vector<point_t>& cooOld);

    std::vector<std::vector<int>>& getRecvDispls()
    {
        return _recvDispls;
    }

    std::vector<std::vector<int>>& getSendDispls()
    {
        return _sendDispls;
    }

    std::vector<MPI_Datatype>& getRecvType()
    {
        return _recvType;
    }

    std::vector<MPI_Datatype>& getSendType()
    {
        return _sendType;
    }

private:
    MPI_Comm _comm;
    MPI_Comm _subComm;
    std::vector<std::vector<int>> _sendDispls;
    std::vector<std::vector<int>> _recvDispls;
    std::vector<MPI_Datatype> _sendType();
    std::vector<MPI_Datatype> _recvType();
};

/// @brief Build send or recv displacement
/// @param absCooNew ending coordinates
/// @param absCooOld starting coordinates
/// @param comm mpi communicator spanning all involved ranks
/// @return send/recv displacement
std::vector<std::vector<int>> getDispls(std::vector<point_t>& absCooNew,
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
/// @param sendDispls
/// @param recvDispls
/// @param redSet
void redDisplsFromSet(std::vector<std::vector<int>>& sendDispls,
   std::vector<std::vector<int>>& recvDispls, std::vector<int>& redSet);

/// @brief Get sub communicator from the reduced set
/// @param redSet
/// @param comm
/// @return
MPI_Comm getSubComm(const std::vector<int>& redSet,
   const MPI_Comm comm = MPI_COMM_WORLD);

/// @brief Get count vector based on displacements
/// @param displs
/// @return
std::vector<int> getCount(const std::vector<std::vector<int>>& displs);


} // namespace Transpose
} // namespace QuICC
