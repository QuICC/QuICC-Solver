/**
 * @file MpiUtils.hpp
 * @brief Methods for mpi enabled transform
 */
#pragma once

// External includes
//
#include <array>
#include <cassert>
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
/// @param sendDispls
/// @param recvDispls
/// @param redSet
void redDisplsFromSet(std::vector<std::vector<int>>& sendDispls,
   std::vector<std::vector<int>>& recvDispls, const std::vector<int>& redSet);

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


/// @brief Container for Mpi communicator and types
/// @tparam TAG
template <class TAG> class Comm
{
public:
   Comm(MPI_Comm comm = MPI_COMM_WORLD) : _comm(comm){};

   /// @brief release resources
   /// like Types and Displacements
   ~Comm()
   {
      for (int r = 0; r < sendType.size(); ++r)
      {
         MPI_Type_free(&sendType[r]);
         MPI_Type_free(&recvType[r]);
      }
   }

   void setComm(const std::vector<point_t>& cooNew,
      const std::vector<point_t>& cooOld)
   {
      _sendDispls = getDispls(cooNew, cooOld);
      _recvDispls = getDispls(cooOld, cooNew);
      auto redSet = getReducedRanksSet(_sendDispls, _recvDispls);
      redDisplsFromSet(_sendDispls, _recvDispls, redSet);
      _subComm = getSubComm(redSet);
   }

   std::vector<std::vector<int>>& getRecvDispls()
   {
      assert(_recvDispls.size() > 0);
      return _recvDispls;
   }

   std::vector<std::vector<int>>& getSendDispls()
   {
      assert(_sendDispls.size() > 0);
      return _sendDispls;
   }

   std::vector<MPI_Datatype>& getRecvType()
   {
      assert(_recvType.size() > 0);
      return _recvType;
   }

   std::vector<MPI_Datatype>& getSendType()
   {
      assert(_sendType.size() > 0);
      return _sendType;
   }

private:
   MPI_Comm _comm;
   MPI_Comm _subComm;
   std::vector<std::vector<int>> _sendDispls;
   std::vector<std::vector<int>> _recvDispls;
   std::vector<MPI_Datatype> _sendType;
   std::vector<MPI_Datatype> _recvType;
};


} // namespace Transpose
} // namespace QuICC
