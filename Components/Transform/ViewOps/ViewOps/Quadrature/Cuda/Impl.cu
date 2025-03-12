
#include <complex>
#include <iostream>
#include <cuda/std/complex>

#include "Impl.hpp"
#include "View/View.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/Quadrature/ViewBatchedMatmulUtils.hpp"
#include "ViewOps/Blas/Cuda/Gemm.hpp"
#include "ViewOps/Blas/Cuda/BlockCGemm.hpp"
#include "ViewOps/Quadrature/Tags.hpp"
#include "ViewOps/ALegendre/Types.hpp"
#include "ViewOps/Worland/Types.hpp"
#include "Cuda/CudaUtil.hpp"

// #define QUICC_USE_NAIVE_CUDA_BATCHED_MATMUL
// #define QUICC_USE_BLOCKED_CUDA_BATCHED_MATMUL
// #define QUICC_USE_VAD_CUDA_BATCHED_MATMUL_REGA
// #define QUICC_USE_VAD_CUDA_BATCHED_MATMUL_REGB
#define QUICC_USE_IBLOCKED_CUDA_MATMUL_WARP_TILE
// #define QUICC_USE_IBLOCKED_CUDA_MATMUL_BLOCK_TILE

namespace QuICC {
namespace Transform {
namespace Quadrature {
namespace Cuda {

namespace details {
using namespace QuICC::Transform::Quadrature::Cuda;
using namespace QuICC::View;

/// cuda kernel of batched matmul with varying sizes
template <class Tout, class Tin, class Top, std::uint16_t Treatment = 0>
__global__ void batchedMatmulKernel(Tout out, const Tin in, const Top op,
   const ViewBase<std::uint32_t> layerIndex,
   const ViewBase<std::uint32_t> layerWidth,
   const ViewBase<std::uint32_t> offSetA, const ViewBase<std::uint32_t> offSetB,
   const ViewBase<std::uint32_t> offSetC)
{
   using IndexType = typename Tin::IndexType;

   const std::size_t l = blockIdx.z;

   // get dimensions
   auto dims = getMatmulDims(out, in, op, layerWidth[l], layerIndex[l]);
   auto M = dims.M;
   auto K = dims.K;
   auto N = dims.N;

   // mats pointers
   auto* c = reinterpret_cast<cuda::std::complex<double>*>(&(out[offSetC[l]]));
   auto* b = reinterpret_cast<cuda::std::complex<double>*>(&(in[offSetB[l]]));
   double* a = &(op[offSetA[l]]);

   // compute
   using alpha_t = typename std::conditional<Treatment != none_m,
      cuda::std::complex<double>, double>::type;
   alpha_t alpha = 1.0;

   if constexpr (Treatment != none_m)
   {
      alpha =
         cuda::std::complex<double>{0.0, static_cast<double>(layerIndex[l])};
      if constexpr (Treatment == diffPhiInt_m)
      {
         alpha = -alpha;
      }
   }
#ifdef QUICC_USE_NAIVE_CUDA_BATCHED_MATMUL
   QuICC::Blas::Cuda::Naive::matmul<double, alpha_t>(c, a, b, M, K, N, alpha);
#endif
#ifdef QUICC_USE_BLOCKED_CUDA_BATCHED_MATMUL
   QuICC::Blas::Cuda::Blocked::matmul<double, alpha_t, 16>(c, a, b, M, K, N,
      alpha);
#endif
#ifdef QUICC_USE_VAD_CUDA_BATCHED_MATMUL_REGA
   QuICC::Blas::Cuda::VadRegA::matmul<double, alpha_t, 64, 16>(c, a, b, M, K, N,
      alpha);
#endif
#ifdef QUICC_USE_VAD_CUDA_BATCHED_MATMUL_REGB
   QuICC::Blas::Cuda::VadRegB::matmul<double, alpha_t, 64, 8>(c, a, b, M, K, N,
      alpha);
#endif
}

/// cuda kernel of batched matmul with varying sizes
template <class Tout, class Tin, class Top, std::uint16_t Treatment = 0>
__global__ void iBlockedMatmulKernel(Tout out, const Tin in, const Top op,
   const ViewBase<std::uint32_t> layerIndex,
   const ViewBase<std::uint32_t> layerWidth,
   const ViewBase<std::uint32_t> offSetA, const ViewBase<std::uint32_t> offSetB,
   const ViewBase<std::uint32_t> offSetC, const ViewBase<std::uint32_t> rowScan,
   const ViewBase<std::uint32_t> colScan, const ViewBase<std::uint32_t> allScan,
   const ViewBase<std::uint32_t> xGrid, const ViewBase<std::uint32_t> yGrid)
{
   using IndexType = typename Tin::IndexType;

   // matrix block ID
   const std::size_t l =
      QuICC::Blas::Cuda::Helper::binary_search_range(allScan, blockIdx.x);
   // get dimensions
   auto dims = getMatmulDims(out, in, op, layerWidth[l], layerIndex[l]);
   auto M = dims.M;
   auto K = dims.K;
   auto N = dims.N;

   // mats pointers
   auto* c = reinterpret_cast<cuda::std::complex<double>*>(&(out[offSetC[l]]));
   auto* b = reinterpret_cast<cuda::std::complex<double>*>(&(in[offSetB[l]]));
   double* a = &(op[offSetA[l]]);

   // last block row and column
   auto lastBlockRow = rowScan[l + 1] - rowScan[l] - 1;
   auto lastBlockCol = colScan[l + 1] - colScan[l] - 1;
   // 2D block coordinates of each block id
   auto blockRow = xGrid[blockIdx.x];
   auto blockCol = yGrid[blockIdx.x];
   QuICC::Blas::Cuda::Helper::BlockIndex lBlock(blockRow, blockCol,
      lastBlockRow, lastBlockCol);

   // compute
   using alpha_t = typename std::conditional<Treatment != none_m,
      cuda::std::complex<double>, double>::type;
   alpha_t alpha = 1.0;
   // compute block matmul
   if constexpr (Treatment != none_m)
   {
      alpha =
         cuda::std::complex<double>{0.0, static_cast<double>(layerIndex[l])};
      if constexpr (Treatment == diffPhiInt_m)
      {
         alpha = -alpha;
      }
   }
#if defined(QUICC_USE_IBLOCKED_CUDA_MATMUL_WARP_TILE)
   QuICC::Blas::Cuda::WarpTile::iBlockGemmDevice<alpha_t>(c, a, b, M, K, N,
      lBlock, alpha, nullptr);
#elif defined(QUICC_USE_IBLOCKED_CUDA_MATMUL_BLOCK_TILE)
   QuICC::Blas::Cuda::BlockTile::iBlockGemmDevice<alpha_t>(c, a, b, M, K, N,
      lBlock, alpha, nullptr);
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
    assert(QuICC::Cuda::isDeviceMemory(out.data()));
    assert(QuICC::Cuda::isDeviceMemory(in.data()));
    assert(QuICC::Cuda::isDeviceMemory(op.data()));

    // batched matmul out = op*in
    assert(op.dims()[0] == out.dims()[0]);
    assert(op.dims()[1] == in.dims()[0]);
    assert(op.dims()[2] == in.dims()[2]);
    assert(op.dims()[2] == out.dims()[2]);

    using IndexType = typename Tin::IndexType;
    using namespace QuICC::Memory;
    using namespace QuICC::View;

    // setup offsets
    if (_layerIndex.data() == nullptr)
    {
        /// \todo move setup to gpu
        auto modsPointers = getModsPointers(out, in, op);

        // copy back to cpu for preprocessing
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
        _layerIndex = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));
        _layerWidth = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));

        // setup view
        ViewBase<IndexType> vLayerIndex(_layerIndex.data(), _layerIndex.size());
        ViewBase<IndexType> vLayerWidth(_layerWidth.data(), _layerWidth.size());

        // setup converters
        tempOnHostMemorySpace converterLI(vLayerIndex, TransferMode::write | TransferMode::block);
        tempOnHostMemorySpace converterLW(vLayerWidth, TransferMode::write);

        IndexType layCtr = 0;
        for (IndexType k = 0; k < modsPointers.size()-1 ; ++k)
        {
            IndexType nCols = modsPointers[k+1] - modsPointers[k];
            // check if layer is populated
            if (nCols > 0)
            {
                vLayerIndex[layCtr] = k;
                vLayerWidth[layCtr] = nCols;
                ++layCtr;
            }
        }

        // alloc
        _offSetA = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers + 1, _mem.get()));
        _offSetB = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers + 1, _mem.get()));
        _offSetC = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers + 1, _mem.get()));
#if defined (QUICC_USE_IBLOCKED_CUDA_MATMUL_WARP_TILE) || defined (QUICC_USE_IBLOCKED_CUDA_MATMUL_BLOCK_TILE)
        _rowScan = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers + 1, _mem.get()));
        _colScan = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers + 1, _mem.get()));
        _allScan = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers + 1, _mem.get()));
#endif

        // setup views
        ViewBase<IndexType> vOffSetA(_offSetA.data(), _offSetA.size());
        ViewBase<IndexType> vOffSetB(_offSetB.data(), _offSetB.size());
        ViewBase<IndexType> vOffSetC(_offSetC.data(), _offSetC.size());
#if defined (QUICC_USE_IBLOCKED_CUDA_MATMUL_WARP_TILE) || defined (QUICC_USE_IBLOCKED_CUDA_MATMUL_BLOCK_TILE)
        ViewBase<IndexType> vRowScan(_rowScan.data(), _rowScan.size());
        ViewBase<IndexType> vColScan(_colScan.data(), _colScan.size());
        ViewBase<IndexType> vAllScan(_allScan.data(), _allScan.size());
#endif


        // setup converters
        tempOnHostMemorySpace converterOA(vOffSetA, TransferMode::write);
        tempOnHostMemorySpace converterOB(vOffSetB, TransferMode::write);
        tempOnHostMemorySpace converterOC(vOffSetC, TransferMode::write);
#if defined (QUICC_USE_IBLOCKED_CUDA_MATMUL_WARP_TILE) || defined (QUICC_USE_IBLOCKED_CUDA_MATMUL_BLOCK_TILE)
        tempOnHostMemorySpace converterRS(vRowScan, TransferMode::write);
        tempOnHostMemorySpace converterCS(vColScan, TransferMode::write);
        tempOnHostMemorySpace converterAS(vAllScan, TransferMode::write);
#endif

        // exclusive scan offsets
        IndexType M, N, K;
        vOffSetA[0] = 0;
        vOffSetB[0] = 0;
        vOffSetC[0] = 0;
#if defined (QUICC_USE_IBLOCKED_CUDA_MATMUL_WARP_TILE) || defined (QUICC_USE_IBLOCKED_CUDA_MATMUL_BLOCK_TILE)
        vRowScan[0] = 0;
        vColScan[0] = 0;
        vAllScan[0] = 0;
#endif
        for (IndexType h = 0; h < nLayers; ++h)
        {
            // get dimensions
            auto dims = getMatmulDims(out, in, op, vLayerWidth[h], vLayerIndex[h]);
            auto M = dims.M;
            auto K = dims.K;
            auto N = dims.N;

            vOffSetA[h+1] = vOffSetA[h] + M*K;
            vOffSetB[h+1] = vOffSetB[h] + K*N;
            vOffSetC[h+1] = vOffSetC[h] + M*N;
#ifdef QUICC_USE_IBLOCKED_CUDA_MATMUL_WARP_TILE
            auto rowBlocks = M / QuICC::Blas::Cuda::WarpTile::BM + (M % QuICC::Blas::Cuda::WarpTile::BM > 0);
            auto colBlocks = N / QuICC::Blas::Cuda::WarpTile::BN + (N % QuICC::Blas::Cuda::WarpTile::BN > 0);

            vRowScan[h+1] = vRowScan[h] + rowBlocks;
            vColScan[h+1] = vColScan[h] + colBlocks;
            vAllScan[h+1] = vAllScan[h] + rowBlocks * colBlocks;
#elif defined(QUICC_USE_IBLOCKED_CUDA_MATMUL_BLOCK_TILE)
            auto rowBlocks = M / QuICC::Blas::Cuda::BlockTile::BM + (M % QuICC::Blas::Cuda::BlockTile::BM > 0);
            auto colBlocks = N / QuICC::Blas::Cuda::BlockTile::BN + (N % QuICC::Blas::Cuda::BlockTile::BN > 0);

            vRowScan[h+1] = vRowScan[h] + rowBlocks;
            vColScan[h+1] = vColScan[h] + colBlocks;
            vAllScan[h+1] = vAllScan[h] + rowBlocks * colBlocks;
#endif
        }

#if defined (QUICC_USE_IBLOCKED_CUDA_MATMUL_WARP_TILE) || defined (QUICC_USE_IBLOCKED_CUDA_MATMUL_BLOCK_TILE)

        _allScanSum = vAllScan[nLayers];
        // alloc
        _xgrid = std::move(
           QuICC::Memory::MemBlock<IndexType>(_allScanSum, _mem.get()));
        _ygrid = std::move(
              QuICC::Memory::MemBlock<IndexType>(_allScanSum, _mem.get()));
        // setup views
        ViewBase<IndexType> vXGrid(_xgrid.data(), _xgrid.size());
        ViewBase<IndexType> vYGrid(_ygrid.data(), _ygrid.size());
        // setup converters

        tempOnHostMemorySpace converterXG(vXGrid, TransferMode::write);
        tempOnHostMemorySpace converterYG(vYGrid, TransferMode::write);
        // generate block cluster
        QuICC::Blas::Cuda::Helper::generate_block_cluster(vXGrid, vYGrid, vRowScan,
           vColScan, vAllScan, nLayers);
#endif
    }

    /// \todo more balanced load distribution
    IndexType M;
    using opLevelType = typename Top::LevelType;
    // if constexpr (ALegendre::is_projector_v<Top>)
    if constexpr (std::is_same_v<opLevelType, CS1RL3D::level> ||
                  std::is_same_v<opLevelType, CS1RL3DJIK::level>)
    {
        M = out.dims()[0]; // longitudinal points
    }
    // else if constexpr (ALegendre::is_integrator_v<Top>)
    else if constexpr (std::is_same_v<opLevelType, S1CLCSC3D::level> ||
                       std::is_same_v<opLevelType, S1CLCSC3DJIK::level>)
    {
        /// \todo store actual max
        M = out.dims()[0] ; // max harmonic order
    }
    // else if constexpr (Worland::Uniform::is_integrator_v<Top>)
    else if constexpr (std::is_same_v<opLevelType, CSL3D::level> ||
                       std::is_same_v<opLevelType, CSL3DJIK::level>)
        M = op.dims()[0]; // radial points/modes
    else
    {
        static_assert(std::is_same_v<opLevelType, void>,
            "backend for these types is not implemented.");
    }

    const IndexType N = _N;
    const IndexType activeLayers = _layerIndex.size();

    // offsets views
    ViewBase<IndexType> layerIndex(_layerIndex.data(), _layerIndex.size());
    ViewBase<IndexType> layerWidth(_layerWidth.data(), _layerWidth.size());
    ViewBase<IndexType> offSetA(_offSetA.data(), _offSetA.size());
    ViewBase<IndexType> offSetB(_offSetB.data(), _offSetB.size());
    ViewBase<IndexType> offSetC(_offSetC.data(), _offSetC.size());
#if defined (QUICC_USE_IBLOCKED_CUDA_MATMUL_WARP_TILE) || defined (QUICC_USE_IBLOCKED_CUDA_MATMUL_BLOCK_TILE)
    ViewBase<IndexType> rowScan(_rowScan.data(), _rowScan.size());
    ViewBase<IndexType> colScan(_colScan.data(), _colScan.size());
    ViewBase<IndexType> allScan(_allScan.data(), _allScan.size());
    // grid views
    ViewBase<IndexType> xGrid(_xgrid.data(), _xgrid.size());
    ViewBase<IndexType> yGrid(_ygrid.data(), _ygrid.size());

    const IndexType allTotal = _allScanSum;

#ifdef QUICC_USE_IBLOCKED_CUDA_MATMUL_WARP_TILE
    dim3 dimBlock(QuICC::Blas::Cuda::WarpTile::GG);
#elif defined(QUICC_USE_IBLOCKED_CUDA_MATMUL_BLOCK_TILE)
    dim3 dimBlock(QuICC::Blas::Cuda::BlockTile::GG);
#endif
    // launch kernel
    details::iBlockedMatmulKernel<Tout, Tin, Top, Treatment>
       <<<allTotal, dimBlock>>>(out, in, op, layerIndex, layerWidth, offSetA,
          offSetB, offSetC, rowScan, colScan, allScan, xGrid, yGrid);
#else

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
    #ifdef QUICC_USE_VAD_CUDA_BATCHED_MATMUL_REGA
    constexpr unsigned int numThreads = 64;
    constexpr unsigned int tileSizeN = 16; // coarsening factor
    blockSize.x = numThreads;
    blockSize.y = 1;
    blockSize.z = 1;
    numBlocks.x = (M + numThreads - 1) / numThreads;
    numBlocks.y = (N + tileSizeN - 1) / tileSizeN;
    numBlocks.z = activeLayers;
    #endif
    #ifdef QUICC_USE_VAD_CUDA_BATCHED_MATMUL_REGB
    constexpr unsigned int numThreads = 64;
    constexpr unsigned int tileSizeM = 8; // coarsening factor
    blockSize.x = numThreads;
    blockSize.y = 1;
    blockSize.z = 1;
    numBlocks.y = (M + tileSizeM - 1) / tileSizeM;
    numBlocks.x = (N + numThreads - 1) / numThreads;
    numBlocks.z = activeLayers;
    #endif

    details::batchedMatmulKernel<Tout, Tin, Top, Treatment>
        <<<numBlocks, blockSize>>>(out, in, op, layerIndex, layerWidth, offSetA, offSetB, offSetC);
#endif

}

// Explicit instantations for AL

// Projectors with row major data
template class ImplOp<ALegendre::physRM_t, ALegendre::modsRM_t, ALegendre::proj_t>;
template class ImplOp<ALegendre::physRM_t, ALegendre::modsRM_t, ALegendre::proj_t, diffPhiPrj_m>;

// Integrators with row major data
template class ImplOp<ALegendre::modsRM_t, ALegendre::physRM_t, ALegendre::int_t>;
template class ImplOp<ALegendre::modsRM_t, ALegendre::physRM_t, ALegendre::int_t, diffPhiInt_m>;


// Projectors with row major data
template class ImplOp<ALegendre::physRM_t, ALegendre::modsRM_t, ALegendre::projRM_t>;
template class ImplOp<ALegendre::physRM_t, ALegendre::modsRM_t, ALegendre::projRM_t, diffPhiPrj_m>;

// Integrators with row major data
template class ImplOp<ALegendre::modsRM_t, ALegendre::physRM_t, ALegendre::intRM_t>;
template class ImplOp<ALegendre::modsRM_t, ALegendre::physRM_t, ALegendre::intRM_t, diffPhiInt_m>;

// Explicit instantations for JW

// Integrators with row major data
template class ImplOp<Worland::Uniform::modsRM_t, Worland::Uniform::physRM_t, Worland::Uniform::int_t>;

// Integrators with row major data
template class ImplOp<Worland::Uniform::modsRM_t, Worland::Uniform::physRM_t, Worland::Uniform::intRM_t>;


} // namespace Cuda
} // namespace Quadrature
} // namespace Transform
} // namespace QuICC
