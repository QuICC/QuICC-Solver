/**
 * @file CommGrouped.hpp
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
#include "ViewOps/Transpose/Mpi/CommUtils.hpp"
#include "Environment/MpiTypes.hpp"
#include "View/View.hpp"
#include "Memory/Memory.hpp"
#include "Memory/Cpu/NewDelete.hpp"
#include "ViewOps/Transpose/Packing.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "Memory/Cuda/Malloc.hpp"
#include "Cuda/CudaUtil.hpp"
#endif

namespace QuICC {
namespace Transpose {
namespace Mpi {

/// @brief Container for Mpi communicator and types
/// to exchange data with MPI_Alltoallv or MPI_Send/MPI_Recv.
/// @tparam TAG implementation tag
template <class TDATA, class TAG = alltoallv_t> class CommGrouped
{
public:
   /// @brief Constructor
   /// @param comm
   CommGrouped(std::shared_ptr<Memory::memory_resource> mem, MPI_Comm comm = MPI_COMM_WORLD) : _mem(mem), _comm(comm){};

   /// @brief release Mpi resources
   ~CommGrouped() = default;

   /// @brief Set the communicator
   /// @param cooNew destination coordinates
   /// @param cooOld source coordinates
   /// @param groupSize number of grouped variables (Views)
   void setComm(const std::vector<point_t>& cooNew,
      const std::vector<point_t>& cooOld, const std::uint64_t groupSize = 1);

   /// @brief Execute the data exchange (communication)
   /// @param out
   /// @param in
   void exchange(std::vector<TDATA*> out, const std::vector<TDATA*> in) const;

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
   /// @brief Entry i specifies the number of elements to send to rank i.
   std::vector<int> _sendCounts;
   /// @brief Entry j specifies the number of elements to receive from rank j.
   std::vector<int> _recvCounts;
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

   /// @brief Max Number of variables to communicate
   std::int64_t _maxGroupSize;

   // #ifdef QUICC_HAS_CUDA_BACKEND
   /// @brief Send buffer displacement for device side packing
   Memory::MemBlock<int> _sendBufferDisplsDevice;
   /// @brief Recv buffer displacement for device side packing
   Memory::MemBlock<int> _recvBufferDisplsDevice;
   /// @brief Send buffer View displacement for device side packing
   View::ViewBase<int> _sendBufferDisplsView;
   /// @brief Recv buffer View displacement for device side packing
   View::ViewBase<int> _recvBufferDisplsView;
   /// @brief Displacement for device side packing
   Memory::MemBlock<int> _sendDisplsDevice;
   /// @brief Displacement for device side packing
   Memory::MemBlock<int> _recvDisplsDevice;
   /// @brief Displacement for device side packing
   View::View<int, View::dense2DRM> _sendDisplsView;
   /// @brief Displacement for device side packing
   View::View<int, View::dense2DRM> _recvDisplsView;
   /// @brief Entry i specifies the number of elements to send to rank i.
   /// Needed for device side packing
   Memory::MemBlock<int> _sendCountsDevice;
   /// @brief Entry j specifies the number of elements to receive from rank j.
   /// Needed for device side packing
   Memory::MemBlock<int> _recvCountsDevice;
   /// @brief View of _sendCountsDevice
   /// Needed for device side packing
   View::ViewBase<int> _sendCountsView;
   /// @brief View of _recvCountsDevice
   /// Needed for device side packing
   View::ViewBase<int> _recvCountsView;
   // #endif

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
   void unPack(TDATA* out, const View::ViewBase<TDATA> buffer) const;

};




template <class TDATA, class TAG>
void CommGrouped<TDATA, TAG>::setComm(const std::vector<point_t>& cooNew,
   const std::vector<point_t>& cooOld, const std::uint64_t groupSize)
{
   _maxGroupSize = groupSize;

   _sendDispls = getDispls(cooNew, cooOld);
   _recvDispls = getDispls(cooOld, cooNew);
   auto redSet = getReducedRanksSet(_sendDispls, _recvDispls);
   redDisplsFromSet(_sendDispls, _recvDispls, redSet);
   _subComm = QuICC::Transpose::Mpi::getSubComm(redSet);
   assert(_subComm != MPI_COMM_NULL);

   auto subRanks = redSet.size();
   _sendCounts = getCount<TAG>(_sendDispls);
   _recvCounts = getCount<TAG>(_recvDispls);

   // Linearized and padded send/recv displacements
   int sendCountsMax = 0;
   int recvCountsMax = 0;
   MPI_Comm_size(_subComm, &_nSubComm);
   for (int i = 0; i < _nSubComm; ++i)
   {
      sendCountsMax = std::max(sendCountsMax, _sendCounts[i]);
      recvCountsMax = std::max(recvCountsMax, _recvCounts[i]);
   }
   _sendDisplsDevice = std::move(Memory::MemBlock<int>(_nSubComm*sendCountsMax , _mem.get()));
   _recvDisplsDevice = std::move(Memory::MemBlock<int>(_nSubComm*recvCountsMax, _mem.get()));

   std::array<std::uint32_t, 2> sendDim {static_cast<std::uint32_t>(_nSubComm), sendCountsMax};
   _sendDisplsView =  View::View<int, View::dense2DRM>({_sendDisplsDevice.data(), _sendDisplsDevice.size()}, sendDim);
   std::array<std::uint32_t, 2> recvDim {static_cast<std::uint32_t>(_nSubComm), recvCountsMax};
   _recvDisplsView =  View::View<int, View::dense2DRM>({_recvDisplsDevice.data(), _recvDisplsDevice.size()}, recvDim);

   // Linearize
   for (int i = 0; i < _nSubComm; ++i)
   {
      for (int j = 0; j < _sendCounts[i]; ++j)
      {
         _sendDisplsView[i*sendCountsMax+j] = _sendDispls[i][j];
      }
      for (int j = 0; j < _recvCounts[i]; ++j)
      {
         _recvDisplsView[i*recvCountsMax+j] = _recvDispls[i][j];
      }
   }

   // Up to here is the same as if the comm were not grouped

   //
   // Buffers
   //

   // scale counts by group size
   for (std::size_t i = 0; i < _nSubComm; ++i)
   {
      _sendCounts[i] *= _maxGroupSize;
      _recvCounts[i] *= _maxGroupSize;
   }

   // Setup send/recv buffers for alltoallv or send/recv
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


   // Buffer offsets
   _sendBufferDisplsView = View::ViewBase<int>(_sendBufferDispls.data(), _sendBufferDispls.size());
   _recvBufferDisplsView = View::ViewBase<int>(_recvBufferDispls.data(), _recvBufferDispls.size());


   // Send/recv Counts
   _sendCountsView = View::ViewBase<int>(_sendCounts.data(), _sendCounts.size());
   _recvCountsView = View::ViewBase<int>(_recvCounts.data(), _recvCounts.size());

   //
   // End Buffers
   //

   _isSetup = true;
}

template <class TDATA, class TAG>
void CommGrouped<TDATA, TAG>::exchange(std::vector<TDATA*> out,
   const std::vector<TDATA*> in) const
{
   if (_subComm != MPI_COMM_NULL)
   {
      // Pack
      #ifdef QUICC_HAS_CUDA_BACKEND
      if(QuICC::Cuda::isDeviceMemory(in[0]))
      {
         Cuda::pack(_sendBufferView, in, _sendCountsView,
            _sendDisplsView, _sendBufferDisplsView);
      }
      else
      #endif
      {
         Cpu::pack(_sendBufferView, in, _sendCountsView,
            _sendDisplsView, _sendBufferDisplsView);
      }

      // CommGrouped
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
         throw std::logic_error("CommGrouped type not implemented");
      }

      // Unpack
      #ifdef QUICC_HAS_CUDA_BACKEND
      if(QuICC::Cuda::isDeviceMemory(out))
      {
         Cuda::unPack(out, _recvBufferView, _recvCountsView,
            _recvDisplsView, _recvBufferDisplsView);
      }
      else
      #endif
      {
         Cpu::unPack(out, _recvBufferView, _recvCountsView,
            _recvDisplsView, _recvBufferDisplsView);
      }
   }
}


} // namespace Mpi
} // namespace Transpose
} // namespace QuICC
