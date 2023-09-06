
#include <complex>
#include <iostream>
#include <cuda/std/complex>

#include "Impl.hpp"
#include "View/View.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/Blas/Cuda/Gemm.hpp"
#include "ViewOps/ALegendre/Tags.hpp"
#include "ViewOps/ALegendre/Types.hpp"
#include "ViewOps/ALegendre/TypeTraits.hpp"
#include "Cuda/CudaUtil.hpp"
#include "Profiler/Interface.hpp"

// #define QUICC_USE_NAIVE_CUDA_BATCHED_MATMUL
// #define QUICC_USE_BLOCKED_CUDA_BATCHED_MATMUL
#define QUICC_USE_VAD_CUDA_BATCHED_MATMUL

namespace QuICC {
namespace Transform {
namespace ALegendre {
namespace Cuda {

namespace details
{
    using namespace QuICC::Transform::ALegendre::Cuda;

    /// cuda kernel of batched matmul with varying sizes
    template<class Tout, class Tin, class Top, std::uint16_t Treatment = 0>
    __global__ void batchedMatmulKernel(Tout out, const Tin in, const Top op,
        const ViewBase<std::uint32_t> harmOrd, const ViewBase<std::uint32_t> cols,
        const ViewBase<std::uint32_t> offSetA, const ViewBase<std::uint32_t> offSetB,
        const ViewBase<std::uint32_t> offSetC)
    {
        using IndexType = typename Tin::IndexType;

        const std::size_t l = blockIdx.z;

        IndexType M, N, K;

        // get dimensions
        if constexpr (is_projector_v<Top>)
        {
            M = out.dims()[0]; // longitudinal points
            K = in.dims()[0] - harmOrd[l]; // current harmonic order
            N = cols[l]; // number of columns
        }
        else if constexpr (is_integrator_v<Top>)
        {
            M = out.dims()[0] - harmOrd[l]; // current harmonic order
            K = in.dims()[0]; // longitudinal points
            N = cols[l]; // number of columns
        }
        else
        {
            static_assert("batched kernel for these types is not implemented.");
        }

        // mats pointers
        auto* c = reinterpret_cast<cuda::std::complex<double>*>(&(out[offSetC[l]]));
        auto* b = reinterpret_cast<cuda::std::complex<double>*>(&(in[offSetB[l]]));
        double* a = &(op[offSetA[l]]);

        // compute
        using alpha_t = typename std::conditional<Treatment == diffPhi_m,
            cuda::std::complex<double>, double>::type;
        alpha_t alpha = 1.0;

        if constexpr (Treatment == diffPhi_m)
        {
            alpha = cuda::std::complex<double>{0.0, static_cast<double>(harmOrd[l])};
            if constexpr (is_integrator_v<Top>)
            {
                alpha = -alpha;
            }
        }
        #ifdef QUICC_USE_NAIVE_CUDA_BATCHED_MATMUL
        QuICC::Blas::Cuda::Naive::matmul<double, alpha_t>(c, a, b, M, K, N, alpha);
        #endif
        #ifdef QUICC_USE_BLOCKED_CUDA_BATCHED_MATMUL
        QuICC::Blas::Cuda::Blocked::matmul<double, alpha_t, 16>(c, a, b, M, K, N, alpha);
        #endif
        #ifdef QUICC_USE_VAD_CUDA_BATCHED_MATMUL
        QuICC::Blas::Cuda::Vad::matmul<double, alpha_t, 64, 16>(c, a, b, M, K, N, alpha);
        #endif
    }
} // namespace details

template<class Tout, class Tin, class Top, std::uint16_t Treatment>
ImplOp<Tout, Tin, Top, Treatment>::ImplOp(std::shared_ptr<QuICC::Memory::memory_resource> mem) : _mem(mem)
{
}


template<class Tout, class Tin, class Top, std::uint16_t Treatment>
void ImplOp<Tout, Tin, Top, Treatment>::applyImpl(Tout& out, const Tin& in, const Top& op)
{
    Profiler::RegionFixture<4> fix("ImplOp::applyImpl");

    assert(QuICC::Cuda::isDeviceMemory(out.data()));
    assert(QuICC::Cuda::isDeviceMemory(in.data()));

    using IndexType = typename Tin::IndexType;

    // setup offsets
    if (_harmOrd.data() == nullptr)
    {
        /// \todo move setup to gpu

        IndexType L; // harmonic order

        ViewBase<IndexType> modsPointers;
        if constexpr (is_projector_v<Top>)
        {
            modsPointers = in.pointers()[1];
            L = in.dims()[0];
        }
        else if constexpr (is_integrator_v<Top>)
        {
            modsPointers = out.pointers()[1];
            L = out.dims()[0];
        }
        else
        {
            static_assert("backend for these types is not implemented.");
        }

        // copy back to cpu for preprocessing
        using namespace QuICC::Memory;
        tempOnHostMemorySpace converterP(modsPointers, TransferMode::read | TransferMode::block);

        _N = 0;
        IndexType nLayers = 0;
        for (IndexType k = 0; k < modsPointers.size()-1 ; ++k)
        {
            IndexType nCols = modsPointers[k+1] - modsPointers[k];
            // check if layer is populated
            if (nCols > 0)
            {
                _N = std::max(_N, nCols);
                ++nLayers;
            }
        }

        // alloc device mem
        _harmOrd = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));
        _cols = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));

        // setup view
        ViewBase<IndexType> vHarmOrd(_harmOrd.data(), _harmOrd.size());
        ViewBase<IndexType> vCols(_cols.data(), _cols.size());

        // setup converters
        tempOnHostMemorySpace converterHO(vHarmOrd, TransferMode::write | TransferMode::block);
        tempOnHostMemorySpace converterC(vCols, TransferMode::write);

        IndexType layCtr = 0;
        for (IndexType k = 0; k < modsPointers.size()-1 ; ++k)
        {
            IndexType nCols = modsPointers[k+1] - modsPointers[k];
            // check if layer is populated
            if (nCols > 0)
            {
                vHarmOrd[layCtr] = k;
                vCols[layCtr] = nCols;
                ++layCtr;
            }
        }

        // alloc
        _offSetA = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));
        _offSetB = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));
        _offSetC = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));

        // setup views
        ViewBase<IndexType> vOffSetA(_offSetA.data(), _offSetA.size());
        ViewBase<IndexType> vOffSetB(_offSetB.data(), _offSetB.size());
        ViewBase<IndexType> vOffSetC(_offSetC.data(), _offSetC.size());

        // setup converters
        tempOnHostMemorySpace converterOA(vOffSetA, TransferMode::write);
        tempOnHostMemorySpace converterOB(vOffSetB, TransferMode::write);
        tempOnHostMemorySpace converterOC(vOffSetC, TransferMode::write);

        // exclusive scan offsets
        IndexType M, N, K;
        vOffSetA[0] = 0;
        vOffSetB[0] = 0;
        vOffSetC[0] = 0;
        for (IndexType h = 0; h < nLayers-1 ; ++h)
        {
            // get dimensions
            if constexpr (is_projector_v<Top>)
            {
                M = out.dims()[0]; // longitudinal points
                K = L - vHarmOrd[h]; // current harmonic order
                N = vCols[h]; // number of columns
            }
            else if constexpr (is_integrator_v<Top>)
            {
                M = L - vHarmOrd[h]; // current harmonic order
                K = in.dims()[0]; // longitudinal points
                N = vCols[h]; // number of columns
            }
            else
            {
                static_assert("backend for these types is not implemented.");
            }

            vOffSetA[h+1] = vOffSetA[h] + M*K;
            vOffSetB[h+1] = vOffSetB[h] + K*N;
            vOffSetC[h+1] = vOffSetC[h] + M*N;
        }
    }


    /// \todo more balanced load distribution
    IndexType M;
    if constexpr (is_projector_v<Top>)
    {
        M = out.dims()[0]; // longitudinal points
    }
    else if constexpr (is_integrator_v<Top>)
    {
        /// \todo store actual max
        M = out.dims()[0] ; // max harmonic order
    }
    else
    {
        static_assert("backend for these types is not implemented.");
    }

    const IndexType N = _N;
    const IndexType activeLayers = _harmOrd.size();

    // offsets views
    ViewBase<IndexType> harmOrd(_harmOrd.data(), _harmOrd.size());
    ViewBase<IndexType> cols(_cols.data(), _cols.size());
    ViewBase<IndexType> offSetA(_offSetA.data(), _offSetA.size());
    ViewBase<IndexType> offSetB(_offSetB.data(), _offSetB.size());
    ViewBase<IndexType> offSetC(_offSetC.data(), _offSetC.size());

    dim3 blockSize;
    dim3 numBlocks;

    #if defined(QUICC_USE_NAIVE_CUDA_BATCHED_MATMUL) || defined(QUICC_USE_BLOCKED_CUDA_BATCHED_MATMUL)
    blockSize.x = 16;
    blockSize.y = 16;
    blockSize.z = 1;
    numBlocks.x = (M + blockSize.x - 1) / blockSize.x;
    numBlocks.y = (N + blockSize.y - 1) / blockSize.y;
    numBlocks.z = activeLayers;
    #endif
    #ifdef QUICC_USE_VAD_CUDA_BATCHED_MATMUL
    constexpr unsigned int numThreads = 64;
    constexpr unsigned int tileSizeN = 16; // coarsening factor
    blockSize.x = numThreads;
    blockSize.y = 1;
    blockSize.z = 1;
    numBlocks.x = (M + blockSize.x - 1) / blockSize.x;
    numBlocks.y = (N + tileSizeN - 1) / tileSizeN;
    numBlocks.z = activeLayers;
    #endif

    details::batchedMatmulKernel<Tout, Tin, Top, Treatment>
        <<<numBlocks, blockSize>>>(out, in, op, harmOrd, cols, offSetA, offSetB, offSetC);

}

// Explicit instantations
using namespace QuICC::Memory;

// Projectors with row major data
template class ImplOp<physRM_t, modsRM_t, proj_t>;
template class ImplOp<physRM_t, modsRM_t, proj_t, diffPhi_m>;

// Integrators with row major data
template class ImplOp<modsRM_t, physRM_t, int_t>;
template class ImplOp<modsRM_t, physRM_t, int_t, diffPhi_m>;


} // namespace Cuda
} // namespace ALegendre
} // namespace Transform
} // namespace QuICC
