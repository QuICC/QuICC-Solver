#include <map>

#include "Graph/Shims/MlirShims.hpp"

namespace QuICC
{
namespace Graph
{

namespace details
{

/// @brief Cache class
/// Computing the size on the device and copying to the host in oreder
/// to allocate the temporary is expensive.
/// We cache the size of the layers in a map based on the metadata pointer.
/// @tparam T type of ptr to cache
template <class T>
class Cache
{
public:
    /// @brief delete copy ctor
    /// @param
    Cache(Cache const&) = delete;
    /// @brief delete assignment op
    /// @param
    void operator=(Cache const&) = delete;
    /// @brief create/return singleton
    /// @return singleton
    static Cache& getInstance();
    /// @brief access the memory resource
    /// @return
    std::map<T*, std::size_t>& getMap();
private:
    /// @brief map pointe to size
    std::map<T*, std::size_t> _ptr2size;
    /// @brief Default ctor
    Cache() = default;
    /// @brief dtor
    ~Cache() = default;
};

template <class T>
Cache<T>& Cache<T>::getInstance()
{
    static Cache instance;
    return instance;
}

template <class T>
std::map<T*, std::size_t>& Cache<T>::getMap()
{
    return _ptr2size;
}

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
    auto& cache = Cache<const std::uint32_t>::getInstance();
    auto it = cache.getMap().find(ptr);
    if (it != cache.getMap().end()) {
        return it->second;
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

    std::uint32_t cumSliceSize = 0;
    std::uint32_t* pCumSliceSize;
    cudaErrChk(cudaMalloc(reinterpret_cast<void**>(&pCumSliceSize), sizeof(std::uint32_t)));
    kernelGetSizeS1CLCSC3DJIK<<<numBlocks, blockSize>>>(pCumSliceSize, ptr, size, lds);
    cudaErrChk(cudaMemcpy(&cumSliceSize, pCumSliceSize, sizeof(std::uint32_t), cudaMemcpyDeviceToHost));
    cudaErrChk(cudaFree(pCumSliceSize));

    // cache the size
    cache.getMap().insert(std::make_pair(ptr, cumSliceSize));
    return cumSliceSize;
}


} // namespace details
} // namespace Graph
} // namespace QuICC