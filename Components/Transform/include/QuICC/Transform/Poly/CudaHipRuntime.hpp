/*
 * @file CudaRuntime.hpp
 * @brief Cuda hip runtime translation layer
 */
#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_CUDAHIPRUNTIME_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_CUDAHIPRUNTIME_HPP


#if defined(__HIP__)
/* #include <hip/hip_runtime.h> */
#ifdef NDEBUG
#undef assert
#define assert(x) ((void)0)
#endif
#define cudaDeviceProp hipDeviceProp_t
#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaErrorInvalidValue hipErrorInvalidValue
#define cudaError_t hipError_t
#define cudaEventCreate hipEventCreate
#define cudaEventDestroy hipEventDestroy
#define cudaEventElapsedTime hipEventElapsedTime
#define cudaEventRecord hipEventRecord
#define cudaEventSynchronize hipEventSynchronize
#define cudaEvent_t hipEvent_t
#define cudaFree hipFree
#define cudaFreeHost hipFreeHost
#define cudaGetDevice hipGetDevice
#define cudaGetDeviceCount hipGetDeviceCount
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaGetErrorName hipGetErrorName
#define cudaGetErrorString hipGetErrorString
#define cudaGetLastError hipGetLastError
#define cudaMalloc hipMalloc
#define cudaMallocHost hipMallocHost
#define cudaMallocManaged hipMallocManaged
#define cudaMemAttachGlobal hipMemAttachGlobal
#define cudaMemcpy hipMemcpy
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemoryTypeDevice hipMemoryTypeDevice
#define cudaPointerAttributes hipPointerAttribute_t
#define cudaPointerGetAttributes hipPointerGetAttributes
#define cudaSetDevice hipSetDevice
#define cudaStreamCreate hipStreamCreate
#define cudaStreamDestroy hipStreamDestroy
#define cudaStreamSynchronize hipStreamSynchronize
#define cudaStream_t hipStream_t
#define cudaSuccess hipSuccess
#define CudaIALegendreOperatorTypes HipIALegendreOperatorTypes
#define CudaIWorlandOperatorTypes HipIWorlandOperatorTypes
#else
/* #include <cuda_runtime.h> */
#endif

#endif //
