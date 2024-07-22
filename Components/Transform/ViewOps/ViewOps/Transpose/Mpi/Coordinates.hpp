/**
 * @file Coordinates.hpp
 * @brief Methods to get absolute coordinate indices
 */
#pragma once

// External includes
//
#include <array>
#include <cassert>
#include <mpi.h>
#include <vector>

// Project includes
//
#include "View/View.hpp"


namespace QuICC
{
namespace View
{

using point_t = std::array<int, 3>;

template <class Tv, class Perm>
std::vector<point_t> getCoo(const Tv& view)
{
    using namespace QuICC::Transpose;
    if constexpr(std::is_same_v<typename Tv::AttributesType, DCCSC3D> && std::is_same_v<Perm, p012_t>)
    {
        std::vector<point_t> coo(view.size());
        auto& pointers = view.pointers()[1];
        auto& indices = view.indices()[1];
        std::size_t itCoo = 0;
        for (std::size_t ptr = 0; ptr < pointers.size()-1; ++ptr)
        {
            for (std::size_t idx = pointers[ptr]; idx < pointers[ptr+1]; ++idx)
            {
                for (std::size_t i = 0; i < view.lds(); ++i)
                {
                    coo[itCoo++] = {
                        static_cast<int>(i),
                        static_cast<int>(indices[idx]),
                        static_cast<int>(ptr)
                        };
                }
            }
        }
        return coo;
    }
    else if constexpr(std::is_same_v<typename Tv::AttributesType, DCCSC3D> && std::is_same_v<Perm, p201_t>)
    {
        std::vector<point_t> coo(view.size());
        auto& pointers = view.pointers()[1];
        auto& indices = view.indices()[1];
        std::size_t itCoo = 0;
        for (std::size_t ptr = 0; ptr < pointers.size()-1; ++ptr)
        {
            for (std::size_t idx = pointers[ptr]; idx < pointers[ptr+1]; ++idx)
            {
                for (std::size_t i = 0; i < view.lds(); ++i)
                {
                    // 2 0 1 -> k i j
                    coo[itCoo++] = {
                        static_cast<int>(ptr),
                        static_cast<int>(i),
                        static_cast<int>(indices[idx])
                        };
                }
            }
        }
        return coo;
    }
    else if constexpr(std::is_same_v<typename Tv::AttributesType, S1CLCSC3D> && std::is_same_v<Perm, p012_t>)
    {
        std::vector<point_t> coo(view.size());
        auto& pointers = view.pointers()[1];
        auto& indices = view.indices()[1];
        std::size_t itCoo = 0;
        for (std::size_t ptr = 0; ptr < pointers.size()-1; ++ptr)
        {
            for (std::size_t idx = pointers[ptr]; idx < pointers[ptr+1]; ++idx)
            {
                std::size_t heightCol = view.dims()[0] - ptr;
                for (std::size_t i = 0; i < heightCol; ++i)
                {
                    coo[itCoo++] = {
                        static_cast<int>(i),
                        static_cast<int>(indices[idx]),
                        static_cast<int>(ptr)
                        };
                }
            }
        }
        return coo;
    }
    else
    {
        throw std::logic_error("getCoo not implemented for this type");
    }
    return {};
}

} // namespace View
} // namespace QuICC

