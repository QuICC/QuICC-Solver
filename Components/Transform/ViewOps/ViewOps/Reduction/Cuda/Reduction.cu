#include <cassert>

#include "Cuda/CudaUtil.hpp"
#include "Profiler/Interface.hpp"
#include "Reduction.hpp"
#include "View/View.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"

namespace QuICC {
namespace Reduction {
namespace Cuda {

using namespace QuICC::View;

namespace details {
// naive implementation of mat reduction
// each thread performs the reduction on a column in the slice
template <class T>
inline __device__ void matReduction(T* out, const T* in, const std::size_t M,
   const std::size_t N)
{
   const std::size_t n = blockIdx.y * blockDim.y + threadIdx.y;

   if (n < N)
   {
      T acc = 0;
      for (std::size_t m = 0; m < M; ++m)
      {
         // in is row major
         auto mn = n + m * N;
         acc += in[mn];
      }
      out[n] = acc;
   }
}

// batched wrapper
// for each z block, extract a slice and call 2D reduction kernel
template <class Tout, class Tin, std::uint32_t Dir>
__global__ void batchedMatReduction(Tout out, const Tin in,
   const ViewBase<std::uint32_t> layerIndex,
   const ViewBase<std::uint32_t> layerWidth,
   const ViewBase<std::uint32_t> offSetOut,
   const ViewBase<std::uint32_t> offSetIn)
{
   static_assert(Dir == 0u &&
                    std::is_same_v<typename Tin::AttributesType, DCCSC3DJIK> &&
                    std::is_same_v<typename Tout::AttributesType, CSC>,
      "Not implemented for other types");

   using IndexType = typename Tin::IndexType;

   const std::size_t l = blockIdx.z;

   // get dimensions
   auto M = in.dims()[0];
   auto N = layerWidth[l];

   // mats pointers
   auto* inPtr = &(in[offSetIn[l]]);
   auto* outPtr = &(out[offSetOut[l]]);

   // compute
   matReduction(outPtr, inPtr, M, N);
}

__device__ void atomicMinDouble(double* addr, double value)
{
    unsigned long long int* addr_as_ull =
        (unsigned long long int*)addr;

    unsigned long long int old = *addr_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(addr_as_ull,
                        assumed,
                        __double_as_longlong(
                            fmin(value, __longlong_as_double(assumed))));
    } while (assumed != old);
}

__device__ inline double warpReduceMin(double val)
{
    for(int offset = warpSize/2; offset > 0; offset /= 2)
        val = fmin(val, __shfl_down_sync(0xffffffff, val, offset));

    return val;
}

__device__ double blockReduceMin(double val)
{
    static __shared__ double shared[4]; // one per warp

    int lane = threadIdx.x % warpSize;
    int wid  = threadIdx.x / warpSize;

    // Step 1: warp reduction
    val = warpReduceMin(val);

    // Step 2: write warp result
    if(lane == 0)
        shared[wid] = val;

    __syncthreads();

    // Step 3: first warp loads warp results
    val = (threadIdx.x < (blockDim.x/warpSize)) ? shared[lane] : 1e308;

    // Step 4: final warp reduction
    if(wid == 0)
        val = warpReduceMin(val);

    return val;
}

template <class Tout, class Tin0, class Tin1, class Tin2, class Tin3, class Tin4, class Tin5>
__global__ void cflMagVel(
    double* r,
    double* dr,
    double* r_ll1,

    const double* vel1,
    const double* vel2,
    const double* vel3,

    const double* mag1,
    const double* mag2,
    const double* mag3,

    int sliceSize,
    int nR,
    int registers_per_thread,

    double mcAlfvenDamping,
    double mcAlfvenScale,

    double* cflRadial,
    double* cflHorizontal)
{
   int i = blockIdx.x * blockDim.x + threadIdx.x;
   int paddedSize = ((sliceSize + blockDim.x*registers_per_thread - 1) / (blockDim.x*registers_per_thread)) * (blockDim.x*registers_per_thread);
   int numBlocksPerSlice = paddedSize / (blockDim.x*registers_per_thread);

   int loc_r = blockIdx.x / (numBlocksPerSlice);
   
    if(loc_r >= nR) return;
    
    double aD;
    double newCfl;

    // Radial CFL
    aD = mcAlfvenDamping / dr[loc_r];
    aD = aD * aD;

    double maxVel = 0.0;

    int localID = (blockIdx.x % numBlocksPerSlice) * (blockDim.x*registers_per_thread) + threadIdx.x;
    for(int j=0;j<registers_per_thread;j++)
    {
        if (localID >= sliceSize) continue;
        double p = mag1[localID + loc_r*sliceSize]*mag1[localID + loc_r*sliceSize] * mcAlfvenScale;
        //if ((localID + loc_r*sliceSize)< 10) printf("gpu %d %e\n", localID, vel1[localID + loc_r*sliceSize]);
        double val =
            p / sqrt(p + aD)
            + abs(vel1[localID + loc_r*sliceSize]);
        //if (val > 1e1) printf("gpu %d %e %e %e\n", localID,  dr[loc_r], val, dr[loc_r] / val);
        if(val > maxVel)
            maxVel = val;

        localID += blockDim.x;
    }
    
    newCfl = dr[loc_r] / maxVel;
    
    newCfl = blockReduceMin(newCfl);
    
    if (threadIdx.x == 0){
        atomicMinDouble(cflRadial, newCfl);

        if(newCfl == cflRadial[0])
            cflRadial[1] = r[loc_r];
    }

    // Horizontal CFL

    aD = mcAlfvenDamping / r_ll1[loc_r];
    aD = aD * aD;

    maxVel = 0.0;

    localID = (blockIdx.x % numBlocksPerSlice) * (blockDim.x*registers_per_thread) + threadIdx.x;

    for(int j=0;j<registers_per_thread;j++)
    {
        if (localID >= sliceSize) continue;

        double p =
            (mag2[localID + loc_r*sliceSize]*mag2[localID + loc_r*sliceSize]
            + mag3[localID + loc_r*sliceSize]*mag3[localID + loc_r*sliceSize])
            * mcAlfvenScale;

        double vel =
            sqrt(
                vel2[localID + loc_r*sliceSize]*vel2[localID + loc_r*sliceSize]
                + vel3[localID + loc_r*sliceSize]*vel3[localID + loc_r*sliceSize]
            );

        double val =
            p / sqrt(p + aD) + vel;
        //if (r_ll1[loc_r] / val < 6e-5) printf("gpu %d %e %e %e\n", localID+ loc_r*sliceSize,  r_ll1[loc_r], val, r_ll1[loc_r] / val);
        if(val > maxVel)
            maxVel = val;

        localID += blockDim.x;
    }

    newCfl = r_ll1[loc_r] / maxVel;
    newCfl = blockReduceMin(newCfl);

    if (threadIdx.x == 0){
        atomicMinDouble(cflHorizontal, newCfl);

        if(newCfl == cflHorizontal[0])
            cflHorizontal[1] = r[loc_r];
    }

}

} // namespace details


template <class Tout, class Tin, std::uint32_t Dir>
Op<Tout, Tin, Dir>::Op(std::shared_ptr<QuICC::Memory::memory_resource> mem) :
    _mem(mem)
{}

template <class Tout, class Tin, std::uint32_t Dir>
void Op<Tout, Tin, Dir>::applyImpl(Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix("Reduction::Cuda::applyImpl");

   assert(QuICC::Cuda::isDeviceMemory(out.data()));
   assert(QuICC::Cuda::isDeviceMemory(in.data()));

   using namespace QuICC::Memory;

   // check types consistency
   static_assert(out.rank() == in.rank() - 1, "input/output rank mismatch");
   static_assert(
      std::is_same_v<typename Tin::ScalarType, typename Tout::ScalarType>,
      "input/output scalar type mismatch");

   if constexpr (Dir == 0u &&
                 std::is_same_v<typename Tin::AttributesType, DCCSC3DJIK> &&
                 std::is_same_v<typename Tout::AttributesType, CSC>)
   {
      // check minimal meta data consistency
      assert(out.pointers()[0].size() == in.pointers()[1].size());
      assert(out.indices()[0].size() == in.indices()[1].size());

      using IndexType = typename Tin::IndexType;
      using namespace QuICC::Memory;

      // setup offsets
      if (_layerIndex.data() == nullptr)
      {
         /// \todo move setup to gpu
         auto pointers = in.pointers()[1];
         assert(pointers.data() != nullptr);

         // copy back to cpu for preprocessing
         tempOnHostMemorySpace converterP(pointers,
            TransferMode::read | TransferMode::block);

         _N = 0;
         IndexType nLayers = 0;
         for (IndexType k = 0; k < pointers.size() - 1; ++k)
         {
            IndexType nCols = pointers[k + 1] - pointers[k];
            // check if layer is populated
            if (nCols > 0)
            {
               _N = std::max(_N, nCols);
               ++nLayers;
            }
         }

         // alloc device mem
         assert(_mem.get() != nullptr);
         _layerIndex =
            std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));
         _layerWidth =
            std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));

         // setup view
         ViewBase<IndexType> vLayerIndex(_layerIndex.data(),
            _layerIndex.size());
         ViewBase<IndexType> vLayerWidth(_layerWidth.data(),
            _layerWidth.size());

         // setup converters
         tempOnHostMemorySpace converterLI(vLayerIndex,
            TransferMode::write | TransferMode::block);
         tempOnHostMemorySpace converterLW(vLayerWidth, TransferMode::write);

         IndexType layCtr = 0;
         for (IndexType k = 0; k < pointers.size() - 1; ++k)
         {
            IndexType nCols = pointers[k + 1] - pointers[k];
            // check if layer is populated
            if (nCols > 0)
            {
               vLayerIndex[layCtr] = k;
               vLayerWidth[layCtr] = nCols;
               ++layCtr;
            }
         }
         // alloc
         _offSetIn =
            std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));
         _offSetOut =
            std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));

         // setup views
         ViewBase<IndexType> vOffSetIn(_offSetIn.data(), _offSetIn.size());
         ViewBase<IndexType> vOffSetOut(_offSetOut.data(), _offSetOut.size());

         // setup converters
         tempOnHostMemorySpace converterOIn(vOffSetIn, TransferMode::write);
         tempOnHostMemorySpace converterOOut(vOffSetOut, TransferMode::write);

         // exclusive scan offsets
         IndexType M, N;
         vOffSetIn[0] = 0;
         vOffSetOut[0] = 0;
         for (IndexType h = 0; h < nLayers - 1; ++h)
         {
            // get dimensions
            auto M = in.dims()[0];
            auto N = vLayerWidth[h];

            vOffSetIn[h + 1] = vOffSetIn[h] + M * N;
            vOffSetOut[h + 1] = vOffSetOut[h] + N;
         }
      }

      /// \todo more balanced load distribution
      IndexType M = in.dims()[0];
      const IndexType N = _N;
      const IndexType activeLayers = _layerIndex.size();

      // offsets views
      ViewBase<IndexType> layerIndex(_layerIndex.data(), _layerIndex.size());
      ViewBase<IndexType> layerWidth(_layerWidth.data(), _layerWidth.size());
      ViewBase<IndexType> offSetIn(_offSetIn.data(), _offSetIn.size());
      ViewBase<IndexType> offSetOut(_offSetOut.data(), _offSetOut.size());

      // setup grid
      dim3 blockSize;
      dim3 numBlocks;

      blockSize.x = 1;
      blockSize.y = 32;
      blockSize.z = 1;
      numBlocks.x = 1;
      numBlocks.y = (N + blockSize.y - 1) / blockSize.y;
      numBlocks.z = activeLayers;

      details::batchedMatReduction<Tout, Tin, Dir><<<numBlocks, blockSize>>>(
         out, in, layerIndex, layerWidth, offSetOut, offSetIn);
   }
   else
   {
      // reduction not implemented
      static_assert(std::is_same_v<typename Tout::AttributesType, void>,
         "Reduction not implemented for this type");
   }
}


template <class Functor, class Tout, class Tin0, class Tin1, class Tin2, class Tin3, class Tin4, class Tin5>
void OpCfl<Functor, Tout, Tin0, Tin1, Tin2, Tin3, Tin4, Tin5>::applyImpl(Tout& outTemp, const Tin0& inVel0, const Tin1& inVel1, const Tin2& inVel2, const Tin3& inMag0, const Tin4& inMag1, const Tin5& inMag2)
{
   Profiler::RegionFixture<4> fix("Reduction::Cuda::applyImpl");

   assert(QuICC::Cuda::isDeviceMemory(inVel0.data()));
   assert(QuICC::Cuda::isDeviceMemory(inMag0.data()));

   using namespace QuICC::Memory;

    // setup grid
    dim3 blockSize;
    dim3 numBlocks;
    int registers_per_thread = 8;

    blockSize.x = 128;
    blockSize.y = 1;
    blockSize.z = 1;
    numBlocks.x = QuICC::Cfl_nR*(((inVel0.dims()[0]*inVel0.dims()[1] + blockSize.x*registers_per_thread - 1)/(blockSize.x*registers_per_thread)));
    numBlocks.y = 1;
    numBlocks.z = 1;

    double radialCfl[2];
    radialCfl[0] = 1e308;
    radialCfl[1] = 0;
    cudaMemcpy(QuICC::newCfl_radial, radialCfl, 2 * sizeof(double),
        cudaMemcpyHostToDevice);

    double horizontalCfl[2];
    horizontalCfl[0] = 1e308;
    horizontalCfl[1] = 0;
    cudaMemcpy(QuICC::newCfl_horizontal, horizontalCfl, 2 * sizeof(double),
        cudaMemcpyHostToDevice);

    details::cflMagVel<Tout, Tin0, Tin1, Tin2, Tin3, Tin4, Tin5><<<numBlocks, blockSize>>>(
        Cfl_r,
        Cfl_dr,
        Cfl_r_ll1,
        inVel0.data(),
        inVel1.data(),
        inVel2.data(),
        inMag0.data(),
        inMag1.data(),
        inMag2.data(),
        inVel0.dims()[0]*inVel0.dims()[1],
        Cfl_nR,
        registers_per_thread,
        Cfl_mcAlfvenDamping,
        Cfl_mcAlfvenScale,
        newCfl_radial,
        newCfl_horizontal);
}

//template <class Tout, class Tin>
//void OpCfl<Tout, Tin>::applyImpl(Tout& out, const Tin& in)
//{
   
//}
// Explicit instantiations
template class Op<View::View<double, CSC>, View::View<double, DCCSC3DJIK>, 0u>;

template class OpCfl<MagVelFunctor<double>, View::View<double, DCCSC3D>, View::View<double, DCCSC3D>, View::View<double, DCCSC3D>, View::View<double, DCCSC3D>, View::View<double, DCCSC3D>, View::View<double, DCCSC3D>, View::View<double, DCCSC3D>>;



} // namespace Cuda
} // namespace Reduction
} // namespace QuICC
