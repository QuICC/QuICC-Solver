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
template <class TDATA, class TAG = void> class Comm
{
public:
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

   void exchange(TDATA* out, const TDATA* in) const
   {
      if (_subComm != MPI_COMM_NULL)
      {
         MPI_Alltoallw(in, _sendCounts.data(), _sDispls.data(),
            _sendType.data(), out, _recvCounts.data(), _rDispls.data(),
            _recvType.data(), _subComm);
      }
   }

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

   bool isSetup() const
   {
      return _isSetup;
   }

private:
   std::vector<std::vector<int>> _sendDispls;
   std::vector<std::vector<int>> _recvDispls;
   std::vector<MPI_Datatype> _sendType;
   std::vector<MPI_Datatype> _recvType;
   std::vector<int> _sendCounts;
   std::vector<int> _recvCounts;
   std::vector<int> _sDispls;
   std::vector<int> _rDispls;
   MPI_Comm _comm;
   MPI_Comm _subComm;
   bool _isSetup = false;
};

} // namespace Mpi
} // namespace Transpose
} // namespace QuICC
