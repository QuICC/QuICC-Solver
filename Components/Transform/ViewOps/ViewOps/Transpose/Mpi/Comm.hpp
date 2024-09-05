/**
 * @file Comm.hpp
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
#include "Environment/MpiTypes.hpp"

namespace QuICC {
namespace Transpose {
namespace Mpi {

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
std::vector<int> getCount(const std::vector<std::vector<int>>& displs);


/// @brief Container for Mpi communicator and types
/// to exchange data with MPI_Alltoallw.
/// @tparam TAG free tag type to implement other
template <class TDATA, class TAG = void> class Comm
{
public:
   /// @brief Constructor
   /// @param comm
   Comm(MPI_Comm comm = MPI_COMM_WORLD) : _comm(comm){};

   /// @brief release Mpi resources
   ~Comm()
   {
      for (std::size_t r = 0; r < _sendType.size(); ++r)
      {
         MPI_Type_free(&_sendType[r]);
         MPI_Type_free(&_recvType[r]);
      }
   }

   /// @brief Execute the data exchange (communication)
   /// @param out
   /// @param in
   void exchange(TDATA* out, const TDATA* in) const
   {
      if (_subComm != MPI_COMM_NULL)
      {
         MPI_Alltoallw(in, _sendCounts.data(), _sDispls.data(),
            _sendType.data(), out, _recvCounts.data(), _rDispls.data(),
            _recvType.data(), _subComm);
      }
   }

   /// @brief Set the communicator
   /// i.e. build the type for the MPI_Alltoallw exchange
   /// @param cooNew destination coordinates
   /// @param cooOld source coordinates
   void setComm(const std::vector<point_t>& cooNew,
      const std::vector<point_t>& cooOld)
   {
      _sendDispls = getDispls(cooNew, cooOld);
      _recvDispls = getDispls(cooOld, cooNew);
      auto redSet = getReducedRanksSet(_sendDispls, _recvDispls);
      redDisplsFromSet(_sendDispls, _recvDispls, redSet);
      _subComm = QuICC::Transpose::Mpi::getSubComm(redSet);
      auto subRanks = redSet.size();
      // Build types
      if (_subComm != MPI_COMM_NULL)
      {
         _sendType.resize(subRanks);
         _recvType.resize(subRanks);

         for (std::size_t r = 0; r < subRanks; ++r)
         {
            MPI_Type_create_indexed_block(_sendDispls[r].size(), 1,
               _sendDispls[r].data(), Environment::MpiTypes::type<TDATA>(),
               &_sendType[r]);
            MPI_Type_commit(&_sendType[r]);
            MPI_Type_create_indexed_block(_recvDispls[r].size(), 1,
               _recvDispls[r].data(), Environment::MpiTypes::type<TDATA>(),
               &_recvType[r]);
            MPI_Type_commit(&_recvType[r]);
         }
      }
      _sendCounts = getCount(_sendDispls);
      _recvCounts = getCount(_recvDispls);
      _sDispls = std::vector<int>(subRanks, 0);
      _rDispls = std::vector<int>(subRanks, 0);
      _isSetup = true;
   }

   /// @brief check if the comm was setup
   /// @return
   bool isSetup() const
   {
      return _isSetup;
   }

private:
   /// @brief Displacement used to create the send types.
   std::vector<std::vector<int>> _sendDispls;
   /// @brief Displacement used to create the recv types.
   std::vector<std::vector<int>> _recvDispls;
   /// @brief Entry i specifies the datatype to use when sending data to rank i.
   std::vector<MPI_Datatype> _sendType;
   /// @brief Entry j specifies the datatype to use when receiving data from
   /// rank j.
   std::vector<MPI_Datatype> _recvType;
   /// @brief Entry i specifies the number of elements to send to rank i.
   std::vector<int> _sendCounts;
   /// @brief Entry j specifies the number of elements to receive from rank j.
   std::vector<int> _recvCounts;
   /// @brief Entry i specifies the displacement (in bytes, offset from sendbuf)
   /// from which to send data to rank i.
   std::vector<int> _sDispls;
   /// @brief Entry j specifies the displacement (in bytes, offset from recvbuf)
   /// to which data from rank j should be written.
   std::vector<int> _rDispls;
   /// @brief All world communicator
   MPI_Comm _comm;
   /// @brief Communicator over which data is to be exchanged.
   MPI_Comm _subComm;
   /// @brief Is comm setup?
   bool _isSetup = false;
};

} // namespace Mpi
} // namespace Transpose
} // namespace QuICC
