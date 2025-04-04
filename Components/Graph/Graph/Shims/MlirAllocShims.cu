#include "Graph/Shims/MlirShims.hpp"
#include "Graph/Shims/CacheLayerSize.hpp"
namespace QuICC
{
namespace Graph
{

namespace details
{

/// @brief kernel to get the size of the layers
/// @param pCumSliceSize
/// @param ptr
/// @param size
/// @param lds
/// @return
__global__ void kernelGetSizeS1CLCSC3DJIK(std::uint32_t* pCumSliceSize, const std::uint32_t* ptr, const std::uint32_t size, const std::uint32_t lds)
{
    /// naive single thread implementation
    std::uint32_t cumSliceSize = 0;
    for (std::uint32_t i = 0; i < size - 1; ++i) {
        auto width = ptr[i+1] - ptr[i];
        auto height = lds - i;
        assert(height > 0);
        cumSliceSize += height * width;
    }
    *pCumSliceSize = cumSliceSize;
}

/// @brief get the cumulative size of the layers
/// @param ptr pointe to ptr metadata
/// @param size
/// @param lds
/// @return size
std::uint32_t getSizeS1CLCSC3DJIK(const std::uint32_t* ptr, const std::uint32_t size, const std::uint32_t lds)
{
    // check if the size was cached
    auto& cache = CacheLayerSize<const std::uint32_t>::getInstance();
    std::uint32_t cumSliceSize = cache.getCachedSizeOrZero(ptr);
    if (cumSliceSize != 0) {
        return cumSliceSize;
    }

    // otherwise, we need to calculate the size

    // setup grid
    dim3 blockSize;
    dim3 numBlocks;

    blockSize.x = 1;
    blockSize.y = 1;
    blockSize.z = 1;
    numBlocks.x = 1;
    numBlocks.y = 1;
    numBlocks.z = 1;

    std::uint32_t* pCumSliceSize;
    cudaErrChk(cudaMalloc(reinterpret_cast<void**>(&pCumSliceSize), sizeof(std::uint32_t)));
    kernelGetSizeS1CLCSC3DJIK<<<numBlocks, blockSize>>>(pCumSliceSize, ptr, size, lds);
    cudaErrChk(cudaMemcpy(&cumSliceSize, pCumSliceSize, sizeof(std::uint32_t), cudaMemcpyDeviceToHost));
    cudaErrChk(cudaFree(pCumSliceSize));

    // cache the size
    cache.addToCache(ptr, cumSliceSize);
    return cumSliceSize;
}


} // namespace details
} // namespace Graph
} // namespace QuICC