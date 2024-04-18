#pragma once


#include <complex>
#include <unordered_map>

#include "ViewOps/ViewMemoryUtils.hpp"

namespace QuICC
{
namespace Graph
{
    template<typename Tdata, typename Tmeta, std::size_t N>
    struct ViewDescriptor {
    Tmeta dims[N];
    Tmeta *pos;
    Tmeta posSize;
    Tmeta *coo;
    Tmeta cooSize;
    Tdata *data;
    Tmeta dataSize;
    };

    using R_DCCSC3D_t = View::View<double, View::DCCSC3D>;
    using C_DCCSC3D_t = QuICC::View::View<std::complex<double>, View::DCCSC3D>;
    using C_S1CLCSC3D_t = QuICC::View::View<std::complex<double>, View::S1CLCSC3D>;

    using view3_t = ViewDescriptor<double, std::uint32_t, 3>;
    using view3_cd_t = ViewDescriptor<std::complex<double>, std::uint32_t, 3>;


    struct ptrAndIdx
    {
        std::vector<std::uint32_t> ptr;
        std::vector<std::uint32_t> idx;
    };

    template <class TagTo, class TagFrom>
    ptrAndIdx denseTransposePtrAndIdx(const std::array<std::uint32_t, 3> dims);

    /// @brief compute pointers and index meta data for fully populated
    /// tensor as ouput of a Transpose operation between AL and JW
    /// @tparam C_DCCSC3D_t
    /// @tparam C_S1CLCSC3D_t
    /// @param dims output dimensions
    /// @return
    template<typename C_DCCSC3D_t, typename C_S1CLCSC3D_t>
    ptrAndIdx denseTransposePtrAndIdx(const std::array<std::uint32_t, 3> dims)
    {
        ptrAndIdx ret;
        // std::uint32_t J = dims[0];
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

} // namespace Graph
} // namespace QuICC


