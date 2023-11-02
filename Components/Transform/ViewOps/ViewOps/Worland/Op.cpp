// External includes
//
#include <iostream>
#include <complex>
#include <type_traits>

// Project includes
//
#include "Op.hpp"
#include "View/View.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/Worland/Impl.hpp"
#include "ViewOps/Worland/Tags.hpp"
#include "ViewOps/Worland/Types.hpp"
#include "ViewOps/Worland/TypeTraits.hpp"
#include "Profiler/Interface.hpp"


namespace QuICC {
namespace Transform {
namespace Worland {

template<class Tout, class Tin, class Top, class Backend>
Op<Tout, Tin, Top, Backend>::Op(span<const typename Top::IndexType> dimensions,
        span<const typename Top::IndexType> layers, std::shared_ptr<QuICC::Memory::memory_resource> mem) :
            _mem(mem)
{
    // forward memory resource if needed
    if constexpr (std::is_constructible_v<Backend>)
    {
        mImpl = std::make_unique<Backend>();
    }
    else
    {
        mImpl = std::make_unique<Backend>(mem);
    }

    auto nLayers = layers.size();

    using IndexType = typename Top::IndexType;
    IndexType M, K, P, pointersSize, indicesSize, metaIdx;
    std::size_t dataSize = 0;

    using LevelType = typename Top::LevelType;
    if constexpr (std::is_same_v<LevelType, CSL3D::level>)
    {
        // uniform truncation projector/integrator
        // dim 0 - Nr - radial points/modes
        // dim 1 - R  - radial modes/points
        // dim 2 - L  - harmonic degree
        M = dimensions[0];
        K = dimensions[1];
        P = dimensions[2];
        pointersSize = 2;
        indicesSize = nLayers;
        metaIdx = 2;
        dataSize = M * K * nLayers;
    }
    else
    {
        throw std::logic_error("ctor for this type is not implemented yet");
    }

    using namespace QuICC::Memory;

    // Alloc op storage
    _opData = MemBlock<typename Top::ScalarType>(dataSize, _mem.get());
    _opPointers = MemBlock<IndexType>(pointersSize, _mem.get());
    _opIndices = MemBlock<IndexType>(indicesSize, _mem.get());

    // Set op view
    ViewBase<IndexType> pointers[_opView.rank()];
    ViewBase<IndexType> indices[_opView.rank()];
    pointers[metaIdx] = ViewBase<IndexType>(_opPointers.data(), _opPointers.size());
    indices[metaIdx] = ViewBase<IndexType>(_opIndices.data(), _opIndices.size());
    _opView = Top(_opData.data(), dataSize, dimensions.data(), pointers, indices);

    // Adapter for device data
    tempOnHostMemorySpace converterP (pointers[metaIdx], TransferMode::write);
    tempOnHostMemorySpace converterI (indices[metaIdx], TransferMode::write);

    // Set up pointers / indices for operator
    if constexpr (std::is_same_v<LevelType, CSL3D::level>)
    {
        // uniform truncation projector/integrator
        pointers[metaIdx][0] = 0;
        pointers[metaIdx][1] = indices[metaIdx].size();

        // Full slice per layer
        for(IndexType p = 0; p < nLayers; ++p)
        {
            // set index
            indices[metaIdx][p] = layers[p];
        }
    }
    else
    {
        throw std::logic_error("ctor for this type is not implemented yet");
    }

}

template<class Tout, class Tin, class Top, class Backend>
void Op<Tout, Tin, Top, Backend>::applyImpl(Tout& out, const Tin& in)
{
    Profiler::RegionFixture<4> fix("Worland::Projector::Op::applyImpl");

    // Apply backend
    mImpl->apply(out, in, _opView);
}

template<class Tout, class Tin, class Top, class Backend>
Top& Op<Tout, Tin, Top, Backend>::getOp()
{
    return _opView;
}

// Explicit instantations

// Projectors and integrators (same data layout) with column major data
template class Op<Uniform::phys_t, Uniform::mods_t, Uniform::projRM_t, Cpu::ImplOp<Uniform::phys_t, Uniform::mods_t, Uniform::projRM_t>>;

// Projectors and integrators (same data layout) with row major data
template class Op<Uniform::physRM_t, Uniform::modsRM_t, Uniform::proj_t, Cpu::ImplOp<Uniform::physRM_t, Uniform::modsRM_t, Uniform::proj_t>>;

// #ifdef QUICC_HAS_CUDA_BACKEND
// template class Op<Uniform::physRM_t, Uniform::modsRM_t, Uniform::proj_t, Cuda::ImplOp<Uniform::physRM_t, Uniform::modsRM_t, Uniform::proj_t>>;
// template class Op<Uniform::physRM_t, Uniform::modsRM_t, Uniform::proj_t, Cuda::ImplOp<Uniform::physRM_t, Uniform::modsRM_t, Uniform::proj_t, diffPhi_m>>;
// #endif



} // namespace Worland
} // namespace Transform
} // namespace QuICC
