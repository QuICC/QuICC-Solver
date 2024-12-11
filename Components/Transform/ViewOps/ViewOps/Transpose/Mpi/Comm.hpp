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
#include <memory>

// Project includes
//
#include "ViewOps/Transpose/Mpi/Tags.hpp"
#include "Environment/MpiTypes.hpp"
#include "View/ViewBase.hpp"
#include "Memory/Memory.hpp"
#include "Memory/Cpu/NewDelete.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "Memory/Cuda/Malloc.hpp"
#include "Cuda/CudaUtil.hpp"
#endif

namespace QuICC {
namespace Transpose {
namespace Mpi {

/// \todo move Mpi utils out of Enviroment and unify
namespace details
{
   inline void mpiAssert(int ierr)
   {
      #ifndef NDEBUG
      if (ierr != MPI_SUCCESS)
      {
         throw  std::runtime_error("Mpi failed.");
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
      if constexpr(std::is_same_v<TAG, alltoallw_t>)
      {
         if (displs[i].size() > 0)
         {
            count[i] = 1;
         }
         else
         {
            count[i] = 0;
         }
      }
      else
      {
         count[i] = displs[i].size();
      }
   }
   return count;
}


/// @brief Container for Mpi communicator and types
/// to exchange data with MPI_Alltoallw, MPI_Alltoallv or MPI_Send/MPI_Recv.
/// @tparam TAG implementation tag
template <class TDATA, class TAG = alltoallv_t> class Comm
{
public:
   /// @brief Constructor
   /// @param comm
   Comm(std::shared_ptr<Memory::memory_resource> mem, MPI_Comm comm = MPI_COMM_WORLD) : _mem(mem), _comm(comm){};

   /// @brief release Mpi resources
   ~Comm();

   /// @brief Set the communicator
   /// i.e. build the type for the MPI_Alltoallw exchange
   /// @param cooNew destination coordinates
   /// @param cooOld source coordinates
   void setComm(const std::vector<point_t>& cooNew,
      const std::vector<point_t>& cooOld);

   /// @brief Execute the data exchange (communication)
   /// @param out
   /// @param in
   void exchange(TDATA* out, const TDATA* in) const;

   /// @brief check if the comm was setup
   /// @return
   bool isSetup() const
   {
      return _isSetup;
   }

private:
   /// @brief Memory resource for buffers
   std::shared_ptr<Memory::memory_resource> _mem;
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
   /// @brief Send buffer for packed comms
   Memory::MemBlock<TDATA> _sendBuffer;
   /// @brief Send buffer view
   View::ViewBase<TDATA> _sendBufferView;
   /// @brief Recv buffer for packed comms
   Memory::MemBlock<TDATA> _recvBuffer;
   /// @brief Recv buffer view
   View::ViewBase<TDATA> _recvBufferView;
   /// @brief Send buffer displacement for packed comms
   std::vector<int> _sendBufferDispls;
   /// @brief Recv buffer displacement for packed comms
   std::vector<int> _recvBufferDispls;
   /// @brief All world communicator
   MPI_Comm _comm;
   /// @brief Communicator over which data is to be exchanged.
   MPI_Comm _subComm;
   /// @brief Sub communicator size
   int _nSubComm;
   /// @brief Is comm setup?
   bool _isSetup = false;

   /// @brief Pack input into buffer for alltoallv and send/recv
   /// @param in
   /// @param buffer
   void pack(View::ViewBase<TDATA> buffer, const TDATA* in) const;

   /// @brief Unpack buffer to output
   /// @param out
   /// @param buffer
   void unPack(TDATA* out, const View::Base<TDATA> buffer) const;

};


template <class TDATA, class TAG>
Comm<TDATA, TAG>::~Comm()
{
   if constexpr(std::is_same_v<TAG, alltoallw_t>)
   {
      for (std::size_t r = 0; r < _sendType.size(); ++r)
      {
         details::mpiAssert(MPI_Type_free(&_sendType[r]));
         details::mpiAssert(MPI_Type_free(&_recvType[r]));
      }
   }
}

template <class TDATA, class TAG>
void Comm<TDATA, TAG>::setComm(const std::vector<point_t>& cooNew,
   const std::vector<point_t>& cooOld)
{
   _sendDispls = getDispls(cooNew, cooOld);
   _recvDispls = getDispls(cooOld, cooNew);
   auto redSet = getReducedRanksSet(_sendDispls, _recvDispls);
   redDisplsFromSet(_sendDispls, _recvDispls, redSet);
   _subComm = QuICC::Transpose::Mpi::getSubComm(redSet);
   auto subRanks = redSet.size();
   _sendCounts = getCount<TAG>(_sendDispls);
   _recvCounts = getCount<TAG>(_recvDispls);
   if (_subComm != MPI_COMM_NULL)
   {
      MPI_Comm_size(_subComm, &_nSubComm);
      if constexpr(std::is_same_v<TAG, alltoallw_t>)
      {
         // Build types for alltoallw
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
         _sDispls = std::vector<int>(subRanks, 0);
         _rDispls = std::vector<int>(subRanks, 0);
      }
      else
      {
         // Setup send/recv buffers for alltoallv or send/recv
         /// \todo gpu buffers
         _sendBufferDispls.resize(_nSubComm+1);
         _recvBufferDispls.resize(_nSubComm+1);
         _sendBufferDispls[0] = 0;
         _recvBufferDispls[0] = 0;
         for (int i = 1; i <= _nSubComm; ++i)
         {
            _sendBufferDispls[i] = _sendBufferDispls[i-1] + _sendCounts[i-1];
            _recvBufferDispls[i] = _recvBufferDispls[i-1] + _recvCounts[i-1];
         }
         _sendBuffer = std::move(Memory::MemBlock<TDATA>(_sendBufferDispls[_nSubComm], _mem.get()));
         _recvBuffer = std::move(Memory::MemBlock<TDATA>(_recvBufferDispls[_nSubComm], _mem.get()));
         // use view so that we have bound checks in debug mode
         _sendBufferView = View::ViewBase<TDATA>(_sendBuffer.data(), _sendBuffer.size());
         _recvBufferView = View::ViewBase<TDATA>(_recvBuffer.data(), _recvBuffer.size());
      }
   }
   _isSetup = true;
}

template <class TDATA, class TAG>
void Comm<TDATA, TAG>::exchange(TDATA* out, const TDATA* in) const
{
   if (_subComm != MPI_COMM_NULL)
   {
      if constexpr(std::is_same_v<TAG, alltoallw_t>)
      {
         details::mpiAssert(MPI_Alltoallw(in, _sendCounts.data(), _sDispls.data(),
            _sendType.data(), out, _recvCounts.data(), _rDispls.data(),
            _recvType.data(), _subComm));
      }
      else
      {
         // Pack
         pack(_sendBufferView, in);

         // Comm
         if constexpr (std::is_same_v<TAG, sendrecv_t>)
         {
            for (int i = 0; i < _nSubComm; ++i)
            {
               details::mpiAssert(MPI_Send(_sendBufferView.data()+_sendBufferDispls[i], _sendCounts[i],
                     Environment::MpiTypes::type<TDATA>(), i, /*tag*/1, _subComm));
            }
            MPI_Status status;
            for (int i = 0; i < _nSubComm; ++i)
            {
               details::mpiAssert(MPI_Recv(_recvBufferView.data()+_recvBufferDispls[i], _recvCounts[i],
                     Environment::MpiTypes::type<TDATA>(), i, /*tag*/1, _subComm, &status));
            }
         }
         else if constexpr (std::is_same_v<TAG, alltoallv_t>)
         {
            details::mpiAssert(MPI_Alltoallv(_sendBufferView.data(), _sendCounts.data(),
               _sendBufferDispls.data(), Environment::MpiTypes::type<TDATA>(),
               _recvBufferView.data(), _recvCounts.data(),
               _recvBufferDispls.data(), Environment::MpiTypes::type<TDATA>(), _subComm));
         }
         else
         {
            throw std::logic_error("Comm type not implemented");
         }

         // Unpack
         unPack(out, _recvBufferView);
      }
   }
}

template <class TDATA, class TAG>
void Comm<TDATA, TAG>::pack(View::ViewBase<TDATA> buffer, const TDATA* in) const
{
   #ifdef QUICC_HAS_CUDA_BACKEND
   assert(Cuda::isDeviceMemory(in) == Cuda::isDeviceMemory(buffer.data()));
   #endif
   for (int i = 0; i < _nSubComm; ++i)
   {
      for (int s = 0; s < _sendCounts[i]; ++s)
      {
         buffer[_sendBufferDispls[i]+s] = *(in + _sendDispls[i][s]);
      }
   }
}

template <class TDATA, class TAG>
void Comm<TDATA, TAG>::unPack(TDATA* out, const View::Base<TDATA> buffer) const
{
   #ifdef QUICC_HAS_CUDA_BACKEND
   assert(Cuda::isDeviceMemory(out) == Cuda::isDeviceMemory(buffer.data()));
   #endif
   for (int i = 0; i < _nSubComm; ++i)
   {
      for (int s = 0; s < _recvCounts[i]; ++s)
      {
         *(out + _recvDispls[i][s]) = buffer[_recvBufferDispls[i]+s];
      }
   }
}

} // namespace Mpi
} // namespace Transpose
} // namespace QuICC
