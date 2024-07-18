/**
 * @file ViewMemoryUtils.hpp
 * @brief Utilities to temporarily transfer memory from the device to host
 * for setup operations that are supported only on the host side
 */
#pragma once

// External includes
//
#include <array>
#include <cstdint>
#include <vector>


// Project includes
//
#include "View/Attributes.hpp"

namespace QuICC {
namespace View {
namespace Index {

    /// @brief logical dimension ordering
    /// these are used to identify that logical data ordering in a view
    enum class dimOrder
    {
        /// @brief xml file dimensions
        v012 = 0,
        /// @brief xml file dimensions
        RThetaPhi = 0,
        /// @brief implicit coriolis spectral spherical ordering
        NLM = 0,
        /// @brief xml file dimensions
        ZXY = 0,

        v021 = 1,
        RPhiTheta = 1,
        /// @brief spectral spherical ordering
        /// input to JW projector, output to AL integrator
        NML = 1,
        ZYX = 1,

        v102 = 2,
        ThetaRPhi = 2,
        /// @brief spectral spherical ordering
        /// input to AL projector, output to FT integrator
        LNM = 2,
        XZY = 2,

        v210 = 3,
        /// @brief physical spherical ordering
        /// input to FT integrator, output to AL projector
        PhiThetaR = 3,
        YXZ = 3,
    };

    struct ptrAndIdx
    {
        std::vector<std::uint32_t> ptr;
        std::vector<std::uint32_t> idx;
    };


    /// @brief Reorder logical dimensions based on order
    /// @param logDims logical dimensions, dimOrder::v012
    /// @param order
    /// @return
    inline std::array<std::uint32_t, 3> getDims(const std::array<std::uint32_t, 3> logDims, const dimOrder order)
    {
        std::array<std::uint32_t, 3> dims;
        switch(order)
        {
            case dimOrder::v012  : dims = {logDims[0], logDims[1], logDims[2]};   break;
            case dimOrder::v021  : dims = {logDims[0], logDims[2], logDims[1]};   break;
            case dimOrder::v102  : dims = {logDims[1], logDims[0], logDims[2]};   break;
            case dimOrder::v210  : dims = {logDims[2], logDims[1], logDims[0]};   break;
            default : throw std::logic_error("case not supported");
        }
        return dims;
    }

    /// @brief compute pointers and index meta data for fully populated
    /// tensor at different stages based on view attributes
    /// @tparam TagViewAttributes
    /// @param logDims logical dimensions, dimOrder::v012
    /// @param order order for which the metadata is created
    /// @return pointer and indices in a struct
    template <class TagViewAttributes>
    ptrAndIdx densePtrAndIdx(const std::array<std::uint32_t, 3> logDims, const dimOrder order = dimOrder::v012)
    {
        assert(false && "not implemented");
        return ptrAndIdx{};
    }


    /// @brief compute pointers and index meta data for fully populated
    /// tensor at different stages based on view attributes
    /// @param logDims logical dimensions, dimOrder::v012
    /// @param order order for which the metadata is created
    /// @return pointer and indices in a struct
    template<>
    ptrAndIdx densePtrAndIdx<DCCSC3D>(const std::array<std::uint32_t, 3> logDims, const dimOrder order)
    {
        std::array<std::uint32_t, 3> dims = getDims(logDims, order);
        auto N = dims[1];
        auto K = dims[2];
        // Populate meta for fully populated tensor
        ptrAndIdx ret;
        ret.ptr.resize(K+1);
        ret.idx.resize(K*N);
        ret.ptr[0] = 0;
        for (std::size_t i = 1; i < ret.ptr.size(); ++i) {
                ret.ptr[i] = ret.ptr[i-1]+N;
        }
        for (std::size_t i = 0; i < ret.idx.size(); ++i) {
                ret.idx[i] = i % N;
        }
        return ret;
    }




} // namespace Index
} // namespace View
} // namespace QuICC