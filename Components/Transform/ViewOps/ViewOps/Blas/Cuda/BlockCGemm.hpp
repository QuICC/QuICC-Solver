/**
 * @file CGemm.hpp
 * @brief Cuda irregular block gemm
 */

#ifndef QUICC_BLAS_BLOCKCGEMM_HPP
#define QUICC_BLAS_BLOCKCGEMM_HPP

// External includes
//
#include <cuda/std/complex>

// Project includes
//
namespace QuICC {

namespace Blas {

namespace Cuda {

namespace Helper {

/* special implementation of the binary search considering found
   if an element is between current and next proc value. */
template <typename V, typename T>
__device__ int binary_search_range(const V scan, const T value)
{
   const int last_index = scan.size() - 1;
   int first_proc = 0;
   int last_proc = last_index;
   auto guess = last_index / 2;

   while (first_proc <= last_proc && first_proc != last_index)
   {
      T current = scan[guess];
      T next = scan[guess + 1];

      if (value >= current && value < next)
      {
         return guess;
      }
      else if (value < current)
      {
         last_proc = guess - 1;
         guess = (first_proc + last_proc + 1) / 2;
      }
      else if (value >= next)
      {
         first_proc = guess + 1;
         guess = (first_proc + last_proc) / 2;
      }
   }

   return -1;
}

struct Abs2Complex
{
   template <typename C> __device__ auto operator()(C v) const
   {
      return (v.real() * v.real()) + (v.imag() * v.imag());
   }
};

template <typename V, typename M, typename I>
void generate_block_cluster(V hxGrid, V hyGrid, M hostRowScan, M hostColScan,
   M hostAllScan, I slowSize)
{
   // Each block is assigned a 2D coordinate  in the matrix
   for (int i = 0; i < slowSize; i++)
   {
      auto rowBlocks = hostRowScan[i + 1] - hostRowScan[i];
      auto colBlocks = hostColScan[i + 1] - hostColScan[i];
      /* Create CUDA Dim Grid (block cluster) as it is not regular.
       * In such way to preserve coalescing by processing column wise*/
      for (int l = 0; l < colBlocks; l++)
      {
         for (int k = 0; k < rowBlocks; k++)
         {
            auto index = hostAllScan[i] + l * rowBlocks + k;
            hxGrid[index] = k;
            hyGrid[index] = l;
         }
      }
   }
}

// block index row and column
template <typename T> struct BlockIndex
{
   T blockRow;
   T blockCol;
   T lastBlockRow;
   T lastBlockCol;

   __device__ BlockIndex(T r, T c, T lbr, T lbc) :
       blockRow(r), blockCol(c), lastBlockRow(lbr), lastBlockCol(lbc)
   {}
};

} // namespace helper

namespace WarpTile {

// warp tile gemm warp and block tile size
constexpr int BM = 32;
constexpr int BN = 32;
constexpr int BK = 16;
constexpr int TM = 4;
constexpr int TN = 2;
constexpr int GG = (BM * BN) / (TM * TN);

// irregular block gemm kernel, default Layout left (L = 0)
template <typename Talpha, typename S, typename D, typename T, typename L,
   typename F = nullptr_t>
__device__ void iBlockGemmDevice(S* Cstart, const D* Astart, const S* Bstart,
   const T M, const T K, const T N, const L lBlock, const Talpha alpha,
   const F& f)
{
   // Thread row and column within Csub
   int row = threadIdx.x / (BN / TN);
   int col = threadIdx.x % (BN / TN);

   auto Axblock = BM;
   auto Ayblock = BN;

   // handle the last row block which can be less than block size.
   if (lBlock.lastBlockRow == lBlock.blockRow)
   {
      Axblock = M - BM * lBlock.lastBlockRow;
   }
   // handle the last col block which can be less than block size.
   if (lBlock.lastBlockCol == lBlock.blockCol)
   {
      Ayblock = N - BN * lBlock.lastBlockCol;
   }

   __shared__ D As[BM * BK];
   __shared__ D Bs1[BK * BN];
   __shared__ D Bs2[BK * BN];

   Astart += lBlock.blockRow * BM * K;
   Bstart += lBlock.blockCol * BN;
   Cstart += lBlock.blockRow * BM * N + lBlock.blockCol * BN;

   const uint innerRowA = threadIdx.x / BK;
   const uint innerColA = threadIdx.x % BK; // warp-level GMEM coalescing
   const uint strideA = GG / BK;
   const uint innerRowB = threadIdx.x / BN;
   const uint innerColB = threadIdx.x % BN; // warp-level GMEM coalescing
   const uint strideB = GG / BN;

   // allocate thread-local cache for results in registerfile
   D threadResults1[TM * TN] = {0.0};
   D threadResults2[TM * TN] = {0.0};
   D regM[TM] = {0.0};
   D reg1N[TN] = {0.0};
   D reg2N[TN] = {0.0};

   auto rs = K % BK;
   auto subMatrixSize = K / BK + (rs > 0);
   // Loop over all the sub-matrices of A and B that are rquired to compute Csub
   // Multiply each pair of sub-matrices together and accumulate the results
   for (int m = 0; m < subMatrixSize; ++m)
   {
      auto ABxyBlock = BK;
      /*handle the last col block which can be less than block size.*/
      if (m == subMatrixSize - 1 && rs > 0)
      {
         ABxyBlock = K - BK * m;
      }

      for (uint loadOffset = 0; loadOffset < BM; loadOffset += strideA)
      {
         /* Handle the last row block */
         if (innerColA < ABxyBlock && innerRowA + loadOffset < Axblock)
         {
            As[(innerRowA + loadOffset) * BK + innerColA] =
               Astart[(innerRowA + loadOffset) * K + innerColA];
         }
      }

      for (uint loadOffset = 0; loadOffset < BK; loadOffset += strideB)
      {
         /* Handle the last col block */
         if (innerRowB + loadOffset < ABxyBlock && innerColB < Ayblock)
         {
            auto bsub = Bstart[(innerRowB + loadOffset) * N + innerColB];
            Bs1[(innerRowB + loadOffset) * BN + innerColB] = bsub.real();
            Bs2[(innerRowB + loadOffset) * BN + innerColB] = bsub.imag();
         }
      }

      /*Synchronize to make sure the sub matrices are loaded
      to shared memory before startingheight the computation*/
      __syncthreads();

      Astart += BK;
      Bstart += BK * N;
      /*Multiply Asub and Bsub together*/
      /* handle the last col and last row block */
      for (int e = 0; e < ABxyBlock; ++e)
      /* for (int e = 0; e < BK; ++e) */
      {
         /* block into registers */
         for (uint i = 0; i < TM; ++i)
         {
            regM[i] = As[(row * TM + i) * BK + e];
         }
         for (uint i = 0; i < TN; ++i)
         {
            reg1N[i] = Bs1[e * BN + col * TN + i];
            reg2N[i] = Bs2[e * BN + col * TN + i];
         }

         for (int resIdxM = 0; resIdxM < TM; ++resIdxM)
         {
            for (int resIdxN = 0; resIdxN < TN; ++resIdxN)
            {
               auto rIdx = resIdxM * TN + resIdxN;
               threadResults1[rIdx] += regM[resIdxM] * reg1N[resIdxN];
               threadResults2[rIdx] += regM[resIdxM] * reg2N[resIdxN];
            }
         }
      }

      /*Synchronize to make sure that the preceding computation is done
      before loading two new sub-matrices of A and B in the next iteration*/
      __syncthreads();
   }

   for (int resIdxM = 0; resIdxM < TM; ++resIdxM)
   {
      for (int resIdxN = 0; resIdxN < TN; ++resIdxN)
      {
         auto rIdx = row * TM + resIdxM;
         auto cIdx = col * TN + resIdxN;
         /* handle the last col and last row block */
         if (rIdx < Axblock && cIdx < Ayblock)
         {
            if constexpr (!std::is_same_v<F, nullptr_t>)
            {
               auto result = S(threadResults1[resIdxM * TN + resIdxN],
                  threadResults2[resIdxM * TN + resIdxN]);
               Cstart[rIdx * N + cIdx] = alpha * f(result);
            }
            else
            {
               auto result = S(threadResults1[resIdxM * TN + resIdxN],
                  threadResults2[resIdxM * TN + resIdxN]);
               Cstart[rIdx * N + cIdx] = alpha * result;
            }
         }
      }
   }
}

} // namespace WarpTile


namespace BlockTile {


// block tile gemm size
constexpr int BLOCK_SIZE = 16;
constexpr int BM = BLOCK_SIZE;
constexpr int BN = BLOCK_SIZE;
constexpr int GG = (BM * BN);

// irregular block gemm kernel, default Layout left (L = 0)
template <typename Talpha, typename S, typename D, typename T, typename L,
   typename F = nullptr_t>
__device__ void iBlockGemmDevice(S* Cstart, const D* Astart, const S* Bstart,
   const T M, const T K, const T N, const L lBlock, const Talpha alpha,
   const F& f)
{
   // Thread row and column within Csub
   int row = threadIdx.x / BLOCK_SIZE;
   int col = threadIdx.x % BLOCK_SIZE;

   auto Axblock = BLOCK_SIZE;
   auto Ayblock = BLOCK_SIZE;

   // handle the last row block which can be less than block size.
   if (lBlock.lastBlockRow == lBlock.blockRow)
   {
      Axblock = M - BLOCK_SIZE * lBlock.lastBlockRow;
   }
   // handle the last col block which can be less than block size.
   if (lBlock.lastBlockCol == lBlock.blockCol)
   {
      Ayblock = N - BLOCK_SIZE * lBlock.lastBlockCol;
   }

   Astart += lBlock.blockRow * BLOCK_SIZE * K;
   Bstart += lBlock.blockCol * BLOCK_SIZE;
   Cstart += lBlock.blockRow * BLOCK_SIZE * N + lBlock.blockCol * BLOCK_SIZE;

   D tmp1 = 0.0;
   D tmp2 = 0.0;

   S Cvalue = 0;
   auto rs = K % BLOCK_SIZE;
   auto subMatrixSize = K / BLOCK_SIZE + (rs > 0);

   // Loop over all the sub-matrices of A and B that are rquired to compute Csub
   // Multiply each pair of sub-matrices together and accumulate the results
   for (int m = 0; m < subMatrixSize; ++m)
   {
      auto ABxyBlock = BLOCK_SIZE;
      /*handle the last col block which can be less than block size.*/
      if (m == subMatrixSize - 1 && rs > 0)
      {
         ABxyBlock = K - BLOCK_SIZE * m;
      }

      __shared__ D As[BLOCK_SIZE][BLOCK_SIZE + 1];
      __shared__ D Bs1[BLOCK_SIZE][BLOCK_SIZE + 1];
      __shared__ D Bs2[BLOCK_SIZE][BLOCK_SIZE + 1];

      /* Handle the last row block */
      if (col < ABxyBlock && row < Axblock)
      {
         As[row][col] = Astart[row * K + col];
      }
      /* Handle the last col block */
      if (row < ABxyBlock && col < Ayblock)
      {
         auto bsub = Bstart[row * N + col];
         Bs1[row][col] = bsub.real();
         Bs2[row][col] = bsub.imag();
      }
      /*Synchronize to make sure the sub matrices are loaded
      to shared memory before startingheight the computation*/
      __syncthreads();

      Astart += BLOCK_SIZE;
      Bstart += BLOCK_SIZE * N;
      /*Multiply Asub and Bsub together*/
      /* handle the last col and last row block */
      if (row < Axblock && col < Ayblock)
      {
         for (int e = 0; e < ABxyBlock; ++e)
         {
            Cvalue += As[row][e] * S(Bs1[e][col], Bs2[e][col]);
         }
      }
      /*Synchronize to make sure that the preceding computation is done
      before loading two new sub-matrices of A and B in the next iteration*/
      __syncthreads();
   }

   /* handle the last col and last row block */
   if (row < Axblock && col < Ayblock)
   {
      S result = 0;
      if constexpr (!std::is_same_v<F, nullptr_t>)
      {
         result = alpha * f(Cvalue);
      }
      else
      {
         result = alpha * Cvalue;
      }
      Cstart[row * N + col] = result;
   }
}

} // namespace BlockTile

} // namespace Cuda
} // namespace Blas
} // namespace QuICC
#endif
