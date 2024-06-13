#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
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
    // buffer size
    data->sizes[0] = lds * idx->sizes[0];
    // alloc data
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
    // buffer size
    data->sizes[0] = lds * idx->sizes[0];
    // alloc data
    details::alloc_ptr(&data->aligned, data->sizes[0], idx->aligned);
};
