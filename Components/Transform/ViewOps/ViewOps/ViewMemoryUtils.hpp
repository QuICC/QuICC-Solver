/**
 * @file ViewMemoryUtils.hpp
 * @brief Utilities to temporarily transfer memory from the device to host
 * for setup operations that are supported only on the host side
 */
#pragma once

// External includes
//

// Project includes
//
#include "View/View.hpp"
#include "Memory/Pensieve.hpp"
#include "Memory/Memory.hpp"
#include "Memory/Cpu/NewDelete.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "Memory/Cuda/Malloc.hpp"
#include "Cuda/CudaUtil.hpp"
#endif

namespace QuICC {
namespace Memory {

namespace TransferMode {
//
// Mode masks
//

/// @brief read mode
constexpr std::uint16_t read = 0b1;

/// @brief write mode
constexpr std::uint16_t write = 0b10;

/// @brief block mode
constexpr std::uint16_t block = 0b100;

} // namespace TransferMode

/// @brief Adapter for temporarily moving Gpu data to Cpu
/// @tparam T
template <class Tview>
class tempOnHostMemorySpace
{
public:
    /// @brief Scalar type
    using ScalarType = typename Tview::ScalarType;

    /// @brief delete default ctor
    tempOnHostMemorySpace() = delete;

    /// @brief ctor referencing view to convert
    /// @param view
    tempOnHostMemorySpace(Tview& view, std::uint16_t mode);

    /// @brief dtor
    ~tempOnHostMemorySpace();
private:
    #ifdef QUICC_HAS_CUDA_BACKEND
    /// @brief temp data storage
    Memory::MemBlock<ScalarType> _dataHost;
    /// @brief store original view data
    View::ViewBase<ScalarType> _dataDevice;
    #endif
    /// @brief store ref to view
    View::ViewBase<ScalarType>& _viewRef{};

    const std::uint16_t _mode;
};


template <class Tview>
tempOnHostMemorySpace<Tview>::tempOnHostMemorySpace(Tview& view, std::uint16_t mode) : _viewRef(view), _mode(mode)
{
    using namespace QuICC::View;
    #ifdef QUICC_HAS_CUDA_BACKEND
    if(QuICC::Cuda::isDeviceMemory(view.data()))
    {
        // wrong space, needs transfer
        static_assert(std::is_same_v<Tview, ViewBase<ScalarType>>,
            "auto transfer implemented only for ViewBase");

        // allocate host memory
        auto& mem = Pensieve<Cpu::NewDelete>::getInstance().getMem();
        _dataHost = MemBlock<ScalarType>(view.size(), &mem);
        // store reference to device memory
        _dataDevice = ViewBase<ScalarType>(view.data(), view.size());
        // redirect view
        view = Tview(_dataHost.data(), _dataHost.size());
        if(_mode & TransferMode::read)
        {
            // transfer data from gpu to temp
            cudaErrChk(cudaMemcpyAsync(_dataHost.data(), _dataDevice.data(),
                _dataDevice.size()*sizeof(ScalarType), cudaMemcpyDeviceToHost));
        }
        if(_mode & TransferMode::block)
        {
            /// \todo make stream aware
            cudaDeviceSynchronize();
        }
    }
    #endif
    // else nothing to do
}

template <class Tview>
tempOnHostMemorySpace<Tview>::~tempOnHostMemorySpace()
{
    #ifdef QUICC_HAS_CUDA_BACKEND
    if(_dataDevice.data() != nullptr)
    {
        // restore view
        _viewRef = Tview(_dataDevice.data(), _dataDevice.size());
        if (_mode & TransferMode::write)
        {
            // transfer temp data to gpu
            cudaErrChk(cudaMemcpyAsync(_dataDevice.data(), _dataHost.data(),
                _dataDevice.size()*sizeof(ScalarType), cudaMemcpyHostToDevice));
            if (_mode & TransferMode::block)
            {
                /// \todo make stream aware
                cudaDeviceSynchronize();
            }
        }
    }
    #endif
    // else nothing to do
}


} // namespace Memory
} // namespace QuICC