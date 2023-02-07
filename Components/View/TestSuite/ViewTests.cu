#include <catch2/catch.hpp>
#include <iostream>

#include "Cuda/CudaUtil.hpp"
#include "Memory/Memory.hpp"
#include "Memory/Cuda/Malloc.hpp"
#include "View/View.hpp"
#include "View/ViewUtils.hpp"


using namespace QuICC::View;
using dense1D = DimLevelType<dense_t>;
using dense2D = DimLevelType<dense_t, dense_t>;

__global__
void add1(View<double, Attributes<dense1D>> out, View<double, Attributes<dense1D>> in)
{
    auto n = out.size();
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index < n) out(index) = in(index) + 1;
}


TEST_CASE("ViewOneDimDenseCuda", "[ViewOneDimDenseCuda]")
{

    constexpr size_t S = 5;
    std::array<double, S> in_h = {1,2,3,4,5};
    std::array<double, S> out_h;
    std::array<double, S> ref = {2,3,4,5,6};

    std::array<std::uint32_t, 1> dimensions {S};

    QuICC::Memory::Cuda::Malloc mem_res;
    QuICC::Memory::MemBlock<double> in_d(S, &mem_res);
    QuICC::Memory::MemBlock<double> out_d(S, &mem_res);

    View<double, Attributes<dense1D>> Vin (in_d, dimensions);
    View<double, Attributes<dense1D>> Vout (out_d, dimensions);

    CHECK(Vin.rank() == 1);
    CHECK(Vin.dims()[0] == S);

    // Copy to gpu
    cudaErrChk(cudaMemcpy(in_d.data(), in_h.data(),
        S*sizeof(double), cudaMemcpyHostToDevice));

    // Invoke kernel 
    const unsigned int blockSize = 32;
    const unsigned int numBlocks = (S + blockSize - 1) / blockSize;
    add1<<<numBlocks, blockSize>>>(Vout, Vin);
    
    // copy back and check
    cudaErrChk(cudaMemcpy(out_h.data(), out_d.data(),
        S*sizeof(double), cudaMemcpyDeviceToHost));

    for (std::size_t i = 0; i < S; ++i)
    {
        CHECK(out_h[i] == ref[i]);
    }

    // test utils
    std::cout << memBlockSize(Vin) << '\n';
}


using sparse1D = DimLevelType<compressed_t>;

__global__
void setI(View<double, Attributes<sparse1D>> data)
{
    auto n = data.pointers()[0][1];
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index < n) data.data()[index] = data.indices()[0][index] + 1;
}

TEST_CASE("ViewOneDimSparseCuda", "[ViewOneDimSparseCuda]")
{
    // Dense data
    constexpr size_t SF = 10;
    std::array<double, SF> fullData = {0,0,3,4,5,0,7,0,0,10};

    // Compressed data
    constexpr size_t S = 5;
    std::array<double, S> data_h = {0,0,0,0,0};
    std::array<std::uint32_t, 1> dimensions {SF};
    std::array<std::vector<std::uint32_t>, 1> pointers_h = {{{0,5}}};
    std::array<std::vector<std::uint32_t>, 1> indices_h = {{{2,3,4,6,9}}};

    // Host view
    View<double, Attributes<sparse1D>> cpuView (data_h, dimensions, pointers_h, indices_h);

    // Get offsets of cpu current and packed pointers
    auto off_h = memOffsets(cpuView);
    auto off_d = memBlockOffsets(cpuView);

    CHECK(off_d[0] == 0);
    CHECK(off_d[1] == S*sizeof(double));
    CHECK(off_d[2] == S*sizeof(double) + 2*sizeof(std::uint32_t));

    // Gpu memory block
    QuICC::Memory::Cuda::Malloc mem_res;
    QuICC::Memory::MemBlock<std::byte> allData_d(memBlockSize(cpuView), &mem_res);

    CHECK(QuICC::Cuda::isDeviceMemory(allData_d.data()) == true);
    CHECK(allData_d.size() == S*sizeof(double) + (2+S)*sizeof(std::uint32_t));

    // Set device pointers (pointers and indices within the device block)
    constexpr size_t rank = 1;
    ViewBase<std::uint32_t> pointers_d[rank];
    pointers_d[0] = ViewBase<std::uint32_t>(reinterpret_cast<std::uint32_t*>(allData_d.data()+off_d[1]), 2);
    ViewBase<std::uint32_t> indices_d[rank];
    indices_d[0] = ViewBase<std::uint32_t>(reinterpret_cast<std::uint32_t*>(allData_d.data()+off_d[2]), S);
    
    CHECK(QuICC::Cuda::isDeviceMemory(pointers_d[0].data()) == true);
    CHECK(QuICC::Cuda::isDeviceMemory(indices_d[0].data()) == true);

    // Device View
    View<double, Attributes<sparse1D>> gpuView (reinterpret_cast<double*>(allData_d.data()), S, dimensions.data(), 
        pointers_d, indices_d);

    // Copy to gpu
    cudaErrChk(cudaMemcpy(allData_d.data(), data_h.data(),
        S*sizeof(double), cudaMemcpyHostToDevice));
    cudaErrChk(cudaMemcpy(gpuView.pointers()[0].data(), pointers_h[0].data(),
        2*sizeof(std::uint32_t), cudaMemcpyHostToDevice));
    cudaErrChk(cudaMemcpy(gpuView.indices()[0].data(), indices_h[0].data(),
        S*sizeof(std::uint32_t), cudaMemcpyHostToDevice));

    // Invoke kernel 
    const unsigned int blockSize = 32;
    const unsigned int numBlocks = (S + blockSize - 1) / blockSize;
    setI<<<numBlocks, blockSize>>>(gpuView);
    
    // Copy back and check
    cudaErrChk(cudaMemcpy(data_h.data(), allData_d.data(),
        S*sizeof(double), cudaMemcpyDeviceToHost));

    for (std::size_t i = 0; i < S; ++i)
    {
        CHECK(data_h[i] == fullData[indices_h[0][i]]);
    }
}


__global__
void addColumn(View<double, Attributes<dense2D>> out, View<double, Attributes<dense2D>> in)
{
    auto m = out.dims()[0];
    auto n = out.dims()[1];
    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if(i < m && j < n) out(i, j) = in(i, j) + j;
}

TEST_CASE("ViewTwoDimDenseColMajCuda", "[ViewTwoDimDenseColMajCuda]")
{
    constexpr size_t M = 3;
    constexpr size_t N = 2;
    constexpr size_t S = M*N;

    std::array<double, S> in_h = {1,1,1,1,1,1};
    std::array<double, S> out_h;
    std::array<double, S> ref = {1,1,1,2,2,2};

    std::array<std::uint32_t, 2> dimensions {M, N};

    QuICC::Memory::Cuda::Malloc mem_res;
    QuICC::Memory::MemBlock<double> in_d(S, &mem_res);
    QuICC::Memory::MemBlock<double> out_d(S, &mem_res);

    View<double, Attributes<dense2D>> Vin (in_d, dimensions);
    View<double, Attributes<dense2D>> Vout (out_d, dimensions);

    CHECK(Vin.rank() == 2);
    CHECK(Vin.dims()[0] == M);
    CHECK(Vin.dims()[1] == N);

    // Copy to gpu
    cudaErrChk(cudaMemcpy(in_d.data(), in_h.data(),
        S*sizeof(double), cudaMemcpyHostToDevice));

    // Invoke kernel 
    dim3 blockSize;
    blockSize.x = 8;
    blockSize.y = 8;
    blockSize.z = 1;
    dim3 numBlocks;
    numBlocks.x = (M + blockSize.x - 1) / blockSize.x;
    numBlocks.y = (N + blockSize.y - 1) / blockSize.y;
    numBlocks.z = 1;

    addColumn<<<numBlocks, blockSize>>>(Vout, Vin);
    
    // copy back and check
    cudaErrChk(cudaMemcpy(out_h.data(), out_d.data(),
        S*sizeof(double), cudaMemcpyDeviceToHost));

    for (std::size_t i = 0; i < S; ++i)
    {
        CHECK(out_h[i] == ref[i]);
    }
}

