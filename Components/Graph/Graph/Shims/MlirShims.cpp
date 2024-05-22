#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "Graph/Types.hpp"


namespace QuICC {
namespace Graph {

/// @brief compute pointers and index meta data for fully populated
/// tensor as ouput of a Transpose operation between AL and JW
/// @tparam C_DCCSC3D_t new buffer type
/// @tparam C_S1CLCSC3D_t producer type
/// @param dims output dimensions
/// @return
template<>
ptrAndIdx denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>(
    const std::array<std::uint32_t, 3> dims)
{
    ptrAndIdx ret;
    std::uint32_t K = dims[1];
    std::uint32_t I = dims[2];

    // row width (with jki) K - ...
    std::vector<std::uint32_t> kLess(I);
    kLess[I-1] = K-1;
    for (std::size_t i = I-1; i > 0; --i)
    {
        if (kLess[i] > 0)
        {
            kLess[i-1] = kLess[i] - 1;
        }
        else
        {
            kLess[i-1] = 0;
        }
    }

    std::size_t layWidthCum = 0;
    ret.ptr.resize(I+1);
    ret.ptr[0] = 0;
    for (std::size_t i = 1; i < I+1; ++i)
    {
        std::uint32_t width = K-kLess[i-1];
        ret.ptr[i] = ret.ptr[i-1] + width;
        layWidthCum += width;
    }

    std::size_t layerIdx = 0;
    ret.idx.resize(layWidthCum);
    for (std::size_t l = 0; l < ret.ptr.size() - 1; ++l)
    {
        auto layerSize = ret.ptr[l+1] - ret.ptr[l];
        for (std::size_t i = 0; i < layerSize; ++i)
        {
            ret.idx[layerIdx+i] = i;
        }
        layerIdx += layerSize;
    }
    return ret;
}

template<>
ptrAndIdx denseTransposePtrAndIdx<C_S1CLCSC3D_t, C_DCCSC3D_t>(const std::array<std::uint32_t, 3> dims)
{
    ptrAndIdx ret;
    std::uint32_t J = dims[1];
    std::uint32_t K = dims[2];

    // fully populated, each layer has all columns
    ret.ptr.resize(K+1);
    ret.ptr[0] = 0;
    for (std::size_t k = 1; k < ret.ptr.size(); ++k)
    {
        ret.ptr[k] = ret.ptr[k-1] + J;
    }

    ret.idx.resize(J*K);
    for (std::size_t i = 1; i < ret.idx.size(); ++i)
    {
        ret.idx[i] = i % J;
    }
    return ret;
}

} // namespace
} // namespace
