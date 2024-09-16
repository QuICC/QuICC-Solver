#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "Graph/Types.hpp"


namespace QuICC {
namespace Graph {

template<>
ptrAndIdx denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>(
    const std::array<std::uint32_t, 3> dims)
{
    ptrAndIdx ret;
    std::uint32_t K = dims[1];
    std::uint32_t I = dims[2];

    // row width (with jki)
    std::uint32_t width = 0;
    std::size_t layWidthCum = 0;
    ret.ptr.resize(I+1);
    ret.ptr[0] = 0;
    for (std::size_t i = 1; i < I+1; ++i)
    {
        width = std::min(width + 1, K);
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

template<>
ptrAndIdx denseTransposePtrAndIdx<C_DCCSC3DJIK_t, C_S1CLCSC3DJIK_t>(
    const std::array<std::uint32_t, 3> dims)
{
    return denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>(dims);
}

template<>
ptrAndIdx denseTransposePtrAndIdx<C_S1CLCSC3DJIK_t, C_DCCSC3DJIK_t>(
    const std::array<std::uint32_t, 3> dims)
{
    return denseTransposePtrAndIdx<C_S1CLCSC3D_t, C_DCCSC3D_t>(dims);
}

} // namespace
} // namespace
