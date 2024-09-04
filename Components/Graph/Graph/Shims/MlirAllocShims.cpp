#include <iostream>
#include <complex>
#include <cassert>

#ifdef QUICC_HAS_CUDA_BACKEND
#include <cuda_runtime.h>
#endif

#include "Graph/Shims/MlirShims.hpp"
#include "Graph/Types.hpp"
#include "View/View.hpp"

using namespace QuICC::Graph;

/// @brief C Interface to MLIR for an allocator of a view data buffer
/// @param data memref
/// @param ptr memref
/// @param idx memref
/// @param lds column height (might be padded for FFT)
extern "C" void _ciface_quiccir_alloc_data_f64_i32_i32_DCCSC3D(
    MemRefDescriptor<double, 1>* data,
    const MemRefDescriptor<std::uint32_t, 1>* ptr,
    const MemRefDescriptor<std::uint32_t, 1>* idx,
    const std::uint64_t lds)
{
    assert(data != nullptr);
    assert(ptr != nullptr);
    assert(idx != nullptr);
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_data_f64_i32_i32_DCCSC3D\n";
    #endif
    // Buffer size
    data->sizes[0] = lds * idx->sizes[0];
    // Alloc buffer
    details::alloc_ptr(&data->aligned, data->sizes[0], idx->aligned);
};

/// @brief C Interface to MLIR for an allocator of a view data buffer
/// @param data memref
/// @param ptr memref
/// @param idx memref
/// @param lds column height (might be padded for FFT)
extern "C" void _ciface_quiccir_alloc_data_complexf64_i32_i32_DCCSC3D(
    MemRefDescriptor<std::complex<double>, 1>* data,
    const MemRefDescriptor<std::uint32_t, 1>* ptr,
    const MemRefDescriptor<std::uint32_t, 1>* idx,
    const std::uint64_t lds)
{
    assert(data != nullptr);
    assert(ptr != nullptr);
    assert(idx != nullptr);
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_data_complexf64_i32_i32_DCCSC3D\n";
    #endif
    // Buffer size
    data->sizes[0] = lds * idx->sizes[0];
    // Alloc buffer
    details::alloc_ptr(&data->aligned, data->sizes[0], idx->aligned);
};

/// @brief C Interface to MLIR for an allocator of a view data buffer
/// @param data memref
/// @param ptr memref
/// @param idx memref
/// @param lds column height (might be padded for FFT)
extern "C" void _ciface_quiccir_alloc_data_complexf64_i32_i32_DCCSC3DJIK(
    MemRefDescriptor<std::complex<double>, 1>* data,
    const MemRefDescriptor<std::uint32_t, 1>* ptr,
    const MemRefDescriptor<std::uint32_t, 1>* idx,
    const std::uint64_t lds)
{
    assert(data != nullptr);
    assert(ptr != nullptr);
    assert(idx != nullptr);
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_data_complexf64_i32_i32_DCCSC3DJIK\n";
    #endif
    // Buffer size
    data->sizes[0] = lds * idx->sizes[0];
    // Alloc buffer
    details::alloc_ptr(&data->aligned, data->sizes[0], idx->aligned);
};

/// @brief C Interface to MLIR for an allocator of a view data buffer
/// @param data memref
/// @param ptr ptr/pos memref
/// @param idx idx/coo memref
/// @param lds column height (might be padded for FFT)
extern "C" void _ciface_quiccir_alloc_data_complexf64_i32_i32_S1CLCSC3D(
    MemRefDescriptor<std::complex<double>, 1>* data,
    const MemRefDescriptor<std::uint32_t, 1>* ptr,
    const MemRefDescriptor<std::uint32_t, 1>* idx,
    const std::uint64_t lds)
{
    assert(data != nullptr);
    assert(ptr != nullptr);
    assert(idx != nullptr);
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_data_complexf64_i32_i32_S1CLCSC3D\n";
    #endif
    // Buffer size
    intptr_t cumSliceSize = 0;
    for (intptr_t i = 0; i < ptr->sizes[0] - 1; ++i) {
        auto width = ptr->aligned[i+1] - ptr->aligned[i];
        auto height = lds - i;
        assert(height > 0);
        cumSliceSize += height * width;
    }
    data->sizes[0] = cumSliceSize;
    // Alloc buffer
    details::alloc_ptr(&data->aligned, data->sizes[0], idx->aligned);
};




/// @brief C Interface to MLIR for an allocator of a view data buffer
/// @param data memref
/// @param ptr ptr/pos memref
/// @param idx idx/coo memref
/// @param lds column height (might be padded for FFT)
extern "C" void _ciface_quiccir_alloc_data_complexf64_i32_i32_S1CLCSC3DJIK(
    MemRefDescriptor<std::complex<double>, 1>* data,
    const MemRefDescriptor<std::uint32_t, 1>* ptr,
    const MemRefDescriptor<std::uint32_t, 1>* idx,
    const std::uint64_t lds)
{
    assert(data != nullptr);
    assert(ptr != nullptr);
    assert(idx != nullptr);
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_data_complexf64_i32_i32_S1CLCSC3DJIK\n";
    #endif
    intptr_t cumSliceSize = 0;
    #ifdef QUICC_HAS_CUDA_BACKEND
    /// \todo replace with mlir side computation of buffer size
    if(QuICC::Cuda::isDeviceMemory(ptr->aligned))
    {
        cumSliceSize = details::getSizeS1CLCSC3DJIK(ptr->aligned, ptr->sizes[0], lds);
    }
    else
    {
    #endif

    // Buffer size
    for (intptr_t i = 0; i < ptr->sizes[0] - 1; ++i) {
        auto width = ptr->aligned[i+1] - ptr->aligned[i];
        auto height = lds - i;
        assert(height > 0);
        cumSliceSize += height * width;
    }
    #ifdef QUICC_HAS_CUDA_BACKEND
    }
    #endif
    data->sizes[0] = cumSliceSize;
    // Alloc buffer
    details::alloc_ptr(&data->aligned, data->sizes[0], idx->aligned);
};