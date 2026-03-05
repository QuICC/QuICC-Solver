/**
 * @file CommGrouped.hpp
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
#include "ViewOps/Transpose/Mpi/CommUtils.hpp"
#include "ViewOps/Transpose/Mpi/Tags.hpp"
#include "ViewOps/Transpose/Packing.hpp"
#include "ViewOps/Transpose/StructArray.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "Cuda/CudaUtil.hpp"
#include "Memory/Cuda/Malloc.hpp"
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
   CommGrouped(std::shared_ptr<Memory::memory_resource> mem,
      MPI_Comm comm = MPI_COMM_WORLD) :
       _mem(mem), _comm(comm){};

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
   template <int SIZE>
   void exchange(structArray<TDATA*, SIZE>& out,
      structArray<const TDATA*, SIZE>& in) const;

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
   /// @brief Send buffer View displacement for device side packing
   View::ViewBase<int> _sendBufferDisplsView;
   /// @brief Recv buffer View displacement for device side packing
   View::ViewBase<int> _recvBufferDisplsView;
   /// @brief Displacement for device side packing
   View::View<int, View::dense2DRM> _sendDisplsView;
   /// @brief Displacement for device side packing
   View::View<int, View::dense2DRM> _recvDisplsView;
   /// @brief View of _sendCountsDevice
   View::ViewBase<int> _sendCountsView;
   /// @brief View of _recvCountsDevice
   View::ViewBase<int> _recvCountsView;
   /// @brief Buffer for linearized displacement for device/host side packing
   Memory::MemBlock<int> _sendDisplsLin;
   /// @brief Buffer for linearized displacement for device/host side packing
   Memory::MemBlock<int> _recvDisplsLin;

   /// @brief Max Number of variables to communicate
   std::int64_t _maxGroupSize;

#ifdef QUICC_HAS_CUDA_BACKEND
   /// @brief Send buffer displacement for device side packing
   Memory::MemBlock<int> _sendBufferDisplsDevice;
   /// @brief Recv buffer displacement for device side packing
   Memory::MemBlock<int> _recvBufferDisplsDevice;
   /// @brief Entry i specifies the number of elements to send to rank i.
   /// Needed for device side packing
   Memory::MemBlock<int> _sendCountsDevice;
   /// @brief Entry j specifies the number of elements to receive from rank j.
   /// Needed for device side packing
   Memory::MemBlock<int> _recvCountsDevice;
#endif

   /// @brief All world communicator
   MPI_Comm _comm;
   /// @brief Communicator over which data is to be exchanged.
   MPI_Comm _subComm;
   /// @brief Sub communicator size
   int _nSubComm;
   /// @brief Is comm setup?
   bool _isSetup = false;
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
   if (_subComm == MPI_COMM_NULL)
   {
      _isSetup = true;
      return;
   }

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
   _sendDisplsLin =
      std::move(Memory::MemBlock<int>(_nSubComm * sendCountsMax, _mem.get()));
   _recvDisplsLin =
      std::move(Memory::MemBlock<int>(_nSubComm * recvCountsMax, _mem.get()));

   std::array<std::uint32_t, 2> sendDim{static_cast<std::uint32_t>(_nSubComm),
      static_cast<std::uint32_t>(sendCountsMax)};
   _sendDisplsView = View::View<int, View::dense2DRM>(
      {_sendDisplsLin.data(), _sendDisplsLin.size()}, sendDim);
   std::array<std::uint32_t, 2> recvDim{static_cast<std::uint32_t>(_nSubComm),
      static_cast<std::uint32_t>(recvCountsMax)};
   _recvDisplsView = View::View<int, View::dense2DRM>(
      {_recvDisplsLin.data(), _recvDisplsLin.size()}, recvDim);

   // Helper views to linearize cpu/gpu data
   View::ViewBase<int> _sendDisplsViewLin(_sendDisplsLin.data(),
      _sendDisplsLin.size());
   View::ViewBase<int> _recvDisplsViewLin(_recvDisplsLin.data(),
      _recvDisplsLin.size());

   // Move temporarly to host
   using namespace QuICC::Memory;
   tempOnHostMemorySpace converterS(_sendDisplsViewLin,
      TransferMode::write | TransferMode::block);
   tempOnHostMemorySpace converterR(_recvDisplsViewLin, TransferMode::write);


   // Linearize
   for (int i = 0; i < _nSubComm; ++i)
   {
      for (int j = 0; j < _sendCounts[i]; ++j)
      {
         _sendDisplsViewLin[i * sendCountsMax + j] = _sendDispls[i][j];
      }
      for (int j = 0; j < _recvCounts[i]; ++j)
      {
         _recvDisplsViewLin[i * recvCountsMax + j] = _recvDispls[i][j];
      }
   }

   // Up to here is the same as if the comm were not grouped

   //
   // Buffers
   //

   // scale counts by group size
   for (int i = 0; i < _nSubComm; ++i)
   {
      _sendCounts[i] *= _maxGroupSize;
      _recvCounts[i] *= _maxGroupSize;
   }

   // Setup send/recv buffers for alltoallv or send/recv
   _sendBufferDispls.resize(_nSubComm + 1);
   _recvBufferDispls.resize(_nSubComm + 1);
   _sendBufferDispls[0] = 0;
   _recvBufferDispls[0] = 0;
   for (int i = 1; i <= _nSubComm; ++i)
   {
      _sendBufferDispls[i] = _sendBufferDispls[i - 1] + _sendCounts[i - 1];
      _recvBufferDispls[i] = _recvBufferDispls[i - 1] + _recvCounts[i - 1];
   }
   _sendBuffer = std::move(
      Memory::MemBlock<TDATA>(_sendBufferDispls[_nSubComm], _mem.get()));
   _recvBuffer = std::move(
      Memory::MemBlock<TDATA>(_recvBufferDispls[_nSubComm], _mem.get()));

   // Use view so that we have bound checks in debug mode
   _sendBufferView =
      View::ViewBase<TDATA>(_sendBuffer.data(), _sendBuffer.size());
   _recvBufferView =
      View::ViewBase<TDATA>(_recvBuffer.data(), _recvBuffer.size());

#ifndef QUICC_HAS_CUDA_BACKEND
   // Buffer offsets
   _sendBufferDisplsView =
      View::ViewBase<int>(_sendBufferDispls.data(), _sendBufferDispls.size());
   _recvBufferDisplsView =
      View::ViewBase<int>(_recvBufferDispls.data(), _recvBufferDispls.size());

   // Send/recv Counts
   _sendCountsView =
      View::ViewBase<int>(_sendCounts.data(), _sendCounts.size());
   _recvCountsView =
      View::ViewBase<int>(_recvCounts.data(), _recvCounts.size());
#else
   // Buffer offsets
   _sendBufferDisplsDevice =
      std::move(Memory::MemBlock<int>(_sendBufferDispls.size(), _mem.get()));
   _recvBufferDisplsDevice =
      std::move(Memory::MemBlock<int>(_recvBufferDispls.size(), _mem.get()));
   _sendBufferDisplsView = View::ViewBase<int>(_sendBufferDisplsDevice.data(),
      _sendBufferDisplsDevice.size());
   _recvBufferDisplsView = View::ViewBase<int>(_recvBufferDisplsDevice.data(),
      _recvBufferDisplsDevice.size());
   // Copy to device
   cudaErrChk(
      cudaMemcpy(_sendBufferDisplsDevice.data(), _sendBufferDispls.data(),
         _sendBufferDispls.size() * sizeof(int), cudaMemcpyHostToDevice));
   cudaErrChk(
      cudaMemcpy(_recvBufferDisplsDevice.data(), _recvBufferDispls.data(),
         _recvBufferDispls.size() * sizeof(int), cudaMemcpyHostToDevice));

   // Send/recv Counts
   _sendCountsDevice =
      std::move(Memory::MemBlock<int>(_sendCounts.size(), _mem.get()));
   _recvCountsDevice =
      std::move(Memory::MemBlock<int>(_recvCounts.size(), _mem.get()));
   _sendCountsView =
      View::ViewBase<int>(_sendCountsDevice.data(), _sendCountsDevice.size());
   _recvCountsView =
      View::ViewBase<int>(_recvCountsDevice.data(), _recvCountsDevice.size());
   // Copy to device
   cudaErrChk(cudaMemcpy(_sendCountsDevice.data(), _sendCounts.data(),
      _sendCounts.size() * sizeof(int), cudaMemcpyHostToDevice));
   cudaErrChk(cudaMemcpy(_recvCountsDevice.data(), _recvCounts.data(),
      _recvCounts.size() * sizeof(int), cudaMemcpyHostToDevice));
#endif

   //
   // End Buffers
   //

   _isSetup = true;
}

template <class TDATA, class TAG>
template <int SIZE>
void CommGrouped<TDATA, TAG>::exchange(structArray<TDATA*, SIZE>& out,
   structArray<const TDATA*, SIZE>& in) const
{
   if (_subComm != MPI_COMM_NULL)
   {
      // Pack
#ifdef QUICC_HAS_CUDA_BACKEND
      if (QuICC::Cuda::isDeviceMemory(in[0]))
      {
         Cuda::pack(_sendBufferView, in, _sendCountsView, _sendDisplsView,
            _sendBufferDisplsView, _maxGroupSize);
      }
      else
#endif
      {
         Cpu::pack(_sendBufferView, in, _sendCountsView, _sendDisplsView,
            _sendBufferDisplsView, _maxGroupSize);
      }

      // CommGrouped
      if constexpr (std::is_same_v<TAG, sendrecv_t>)
      {
         for (int i = 0; i < _nSubComm; ++i)
         {
            details::mpiAssert(MPI_Send(
               _sendBufferView.data() + _sendBufferDispls[i], _sendCounts[i],
               Environment::MpiTypes::type<TDATA>(), i, /*tag*/ 1, _subComm));
         }
         MPI_Status status;
         for (int i = 0; i < _nSubComm; ++i)
         {
            details::mpiAssert(
               MPI_Recv(_recvBufferView.data() + _recvBufferDispls[i],
                  _recvCounts[i], Environment::MpiTypes::type<TDATA>(), i,
                  /*tag*/ 1, _subComm, &status));
         }
      }
      else if constexpr (std::is_same_v<TAG, alltoallv_t>)
      {
         /*details::mpiAssert(MPI_Alltoallv(_sendBufferView.data(),
            _sendCounts.data(), _sendBufferDispls.data(),
            Environment::MpiTypes::type<TDATA>(), _recvBufferView.data(),
            _recvCounts.data(), _recvBufferDispls.data(),
            Environment::MpiTypes::type<TDATA>(), _subComm));*/

         std::vector<MPI_Request> reqs;
         reqs.reserve(2 * _nSubComm);
         int rank;
         MPI_Comm_rank(_subComm, &rank);
         int typesize;
        MPI_Type_size(Environment::MpiTypes::type<TDATA>(), &typesize);
         cudaMemcpyAsync(_recvBufferView.data() + _recvBufferDispls[rank],
                   _sendBufferView.data() + _sendBufferDispls[rank],
                   _sendCounts[rank] * typesize, cudaMemcpyDeviceToDevice, 0);

         for (int i = 0; i < _nSubComm; ++i)
         {
            if ((_recvCounts[i] > 0) && (i!=rank))
            {
               MPI_Request r;
               MPI_Irecv(_recvBufferView.data() + _recvBufferDispls[i],
                  _recvCounts[i], Environment::MpiTypes::type<TDATA>(), i, 0,
                  _subComm, &r);
               reqs.push_back(r);
            }
         }
         for (int i = 0; i < _nSubComm; ++i)
         {
            if ((_sendCounts[i] > 0) && (i!=rank))
            {
               MPI_Request r;
               MPI_Isend(_sendBufferView.data() + _sendBufferDispls[i],
                  _sendCounts[i], Environment::MpiTypes::type<TDATA>(), i, 0,
                  _subComm, &r);
               reqs.push_back(r);
            }
         }
         MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

         /* MPI_Request reqs[2];
         int rank;
         MPI_Comm_rank(_subComm, &rank);
         for (int i = 0; i < _nSubComm; ++i)
         {
             if (i == 0)
             {
                cudaMemcpyAsync(_recvBufferView.data() + _recvBufferDispls[rank],
                   _sendBufferView.data() + _sendBufferDispls[rank],
                   _sendCounts[rank] * sizeof(Environment::MpiTypes::type<TDATA>()), cudaMemcpyDeviceToDevice, 0);
             }
             else
             {
                int next = (rank + i) % _nSubComm;
                int prev = (rank - i + _nSubComm) % _nSubComm;
                MPI_Irecv(_recvBufferView.data() + _recvBufferDispls[prev],
                   _recvCounts[prev], Environment::MpiTypes::type<TDATA>(),
                   prev, 0, _subComm, &reqs[0]);
                MPI_Isend(_sendBufferView.data() + _sendBufferDispls[next],
                   _sendCounts[next], Environment::MpiTypes::type<TDATA>(),
                   next, 0, _subComm, &reqs[1]);
                MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);
             }
         }*/
         cudaDeviceSynchronize();
      }
      else
      {
         throw std::logic_error("CommGrouped type not implemented");
      }

// Unpack
#ifdef QUICC_HAS_CUDA_BACKEND
      if (QuICC::Cuda::isDeviceMemory(out[0]))
      {
          //Cuda::unPack(out, _sendBufferView, _sendCountsView, _sendDisplsView,
          //  _sendDisplsView, _maxGroupSize);
         Cuda::unPack(out, _recvBufferView, _recvCountsView, _recvDisplsView,
            _recvBufferDisplsView, _maxGroupSize);
      }
      else
#endif
      {
         Cpu::unPack(out, _recvBufferView, _recvCountsView, _recvDisplsView,
            _recvBufferDisplsView, _maxGroupSize);
      }
   }
}


} // namespace Mpi
} // namespace Transpose
} // namespace QuICC
