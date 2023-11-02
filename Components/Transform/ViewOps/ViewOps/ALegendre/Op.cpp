
#include <iostream>
#include <complex>
#include <type_traits>

#include "Op.hpp"
#include "View/View.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/ALegendre/Impl.hpp"
#include "ViewOps/ALegendre/Tags.hpp"
#include "ViewOps/ALegendre/Types.hpp"
#include "ViewOps/ALegendre/TypeTraits.hpp"
#include "Profiler/Interface.hpp"


namespace QuICC {
namespace Transform {
namespace ALegendre {


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
    if constexpr (is_integrator_v<Top>)
    {
        /// dim 0 - L  - harmonic degree
        /// dim 1 - Nl - longitudinal points
        /// dim 2 - M  - harmonic order
        K = dimensions[0];
        M = dimensions[1];
        P = dimensions[2];
        pointersSize = P+1;
        indicesSize = nLayers*M;
        metaIdx = 1;
    }
    else if constexpr (is_projector_v<Top>)
    {
        /// dim 0 - Nl - longitudinal points
        /// dim 1 - L  - harmonic degree
        /// dim 2 - M  - harmonic order
        M = dimensions[0];
        K = dimensions[1];
        P = dimensions[2];
        pointersSize = 2;
        indicesSize = nLayers;
        metaIdx = 2;
    }
    else
    {
        throw std::logic_error("ctor for this type is not implemented yet");
    }

    std::size_t varSize = 0;
    for(IndexType i = 0; i < nLayers; ++i)
    {
        auto nPoly = K - layers[i];
        varSize += nPoly;
    }
    std::size_t dataSize = M * varSize;

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
    if constexpr (is_integrator_v<Top>)
    {
        pointers[metaIdx][0] = 0;
        IndexType nCols = 0;
        IndexType layerCtr = 0;
        for (IndexType i = 0; i < P; ++i)
        {
            IndexType blockWidth = 0;
            if (static_cast<IndexType>(layers[layerCtr] == i))
            {
                blockWidth = M;
                ++layerCtr;
            }

            for (IndexType idx = 0; idx < blockWidth; idx++)
            {
                indices[metaIdx][nCols+idx] = idx;
            }

            nCols += blockWidth;
            pointers[metaIdx][i+1] = nCols;
        }
    }
    else if constexpr (is_projector_v<Top>)
    {
        pointers[metaIdx][0] = 0;
        pointers[metaIdx][1] = indices[metaIdx].size();

        // Copy operators matrix by matrix to view
        for(IndexType p = 0; p < nLayers; ++p)
        {
            // set index
            indices[metaIdx][p] = layers[p];
        }
    }
}

template<class Tout, class Tin, class Top, class Backend>
void Op<Tout, Tin, Top, Backend>::applyImpl(Tout& out, const Tin& in)
{
    Profiler::RegionFixture<4> fix("ALegendre::Projector::Op::applyImpl");

    // Apply backend
    mImpl->apply(out, in, _opView);
}

template<class Tout, class Tin, class Top, class Backend>
Top& Op<Tout, Tin, Top, Backend>::getOp()
{
    return _opView;
}

// Explicit instantations
using namespace QuICC::Memory;

// Projectors with column major data
template class Op<phys_t, mods_t, projRM_t, Cpu::ImplOp<phys_t, mods_t, projRM_t>>;
template class Op<phys_t, mods_t, projRM_t, Cpu::ImplOp<phys_t, mods_t, projRM_t, diffPhi_m>>;

// Projectors with row major data
template class Op<physRM_t, modsRM_t, proj_t, Cpu::ImplOp<physRM_t, modsRM_t, proj_t>>;
template class Op<physRM_t, modsRM_t, proj_t, Cpu::ImplOp<physRM_t, modsRM_t, proj_t, diffPhi_m>>;

#ifdef QUICC_HAS_CUDA_BACKEND
template class Op<physRM_t, modsRM_t, proj_t, Cuda::ImplOp<physRM_t, modsRM_t, proj_t>>;
template class Op<physRM_t, modsRM_t, proj_t, Cuda::ImplOp<physRM_t, modsRM_t, proj_t, diffPhi_m>>;
#endif

// Integrators with column major data
template class Op<mods_t, phys_t, intRM_t, Cpu::ImplOp<mods_t, phys_t, intRM_t>>;
template class Op<mods_t, phys_t, intRM_t, Cpu::ImplOp<mods_t, phys_t, intRM_t, diffPhi_m>>;

// Integrators with row major data
template class Op<modsRM_t, physRM_t, int_t, Cpu::ImplOp<modsRM_t, physRM_t, int_t>>;
template class Op<modsRM_t, physRM_t, int_t, Cpu::ImplOp<modsRM_t, physRM_t, int_t, diffPhi_m>>;

#ifdef QUICC_HAS_CUDA_BACKEND
template class Op<modsRM_t, physRM_t, int_t, Cuda::ImplOp<modsRM_t, physRM_t, int_t>>;
template class Op<modsRM_t, physRM_t, int_t, Cuda::ImplOp<modsRM_t, physRM_t, int_t, diffPhi_m>>;
#endif

} // namespace ALegendre
} // namespace Transform
} // namespace QuICC
