/**
 * @file PIALegendreOperatorUtils.hpp
 * @brief Associated Legendre based operator cuda irregular block gemm utils
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_CUDAIALEGENDREOPERATORGEMMUTILS_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_CUDAIALEGENDREOPERATORGEMMUTILS_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

#ifdef QUICC_USE_KOKKOS_CUDA
#include "QuICC/Transform/Poly/KokkosUtils.hpp"
#include <cuComplex.h>
#include "CudaIALegendreOperatorTypes.hpp"
#endif

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

#ifdef QUICC_USE_KOKKOS_CUDA
/* using DataType = cuDoubleComplex; */

//GEMM UTILS
//
// Thread block size
#define BLOCK_SIZE 16
/* #define BLOCK_SIZE 32 */

// Get a matrix element, layout left
template <Integer Layout, typename M>
__device__ std::enable_if_t<Layout == 0, typename M::Scalar> GetElement(
   const M A, int row, int col) {
   return A.elements[col * A.stride + row];
}

// Get a matrix element, layout right
template <Integer Layout, typename M>
__device__ std::enable_if_t<Layout == 1, typename M::Scalar> GetElement(
   const M A, int row, int col) {
   return A.elements[row * A.stride + col];
}

// Set a matrix element, layout left, cuda complex
template <typename M>
__device__ void SetElement(M A, int row, int col, cuDoubleComplex value) {
   A.elements[col * A.stride + row].real() = cuCreal(value);
   A.elements[col * A.stride + row].imag() = cuCimag(value);
}

// Set a matrix element
template <Integer Layout, typename M>
__device__ std::enable_if_t<Layout == 0, void> SetElement(
   M A, int row, int col, typename M::Scalar value) {
   A.elements[col * A.stride + row] = value;
}

// Set a matrix element
template <Integer Layout, typename M>
__device__ std::enable_if_t<Layout == 1, void> SetElement(
   M A, int row, int col, typename M::Scalar value) {
   A.elements[row * A.stride + col] = value;
}

template <typename M> __device__ M GetStartMatrix(M A, int offset, int stride, int width) {
   M Asub;
   Asub.width = width;
   Asub.height = stride;
   Asub.block_stride = stride;
   Asub.stride = A.stride;
   Asub.elements = &A.elements[offset];
   return Asub;
}

template <typename M> __device__ M GetStartMatrix(M A, int offset, int stride) {
   M Asub;
   Asub.height = stride;
   Asub.block_stride = stride;
   Asub.stride = A.stride;
   Asub.elements = &A.elements[offset];
   return Asub;
}

// Get the BLOCK_SIZExBLOCK_SIZE sub-matrix Asub of A that is
// located col sub-matrices to the right and row sub-matrices down
// from the upper-left corner of A
template <Integer Layout, typename M>
__device__ std::enable_if_t<Layout == 0, M> GetSubMatrix(
   M A, int row, int col) {
   M Asub;
   Asub.width = BLOCK_SIZE;
   Asub.height = BLOCK_SIZE;
   Asub.stride = A.stride;
   Asub.elements = &A.elements[A.stride * BLOCK_SIZE * col + BLOCK_SIZE * row];
   return Asub;
}

// Get the BLOCK_SIZExBLOCK_SIZE sub-matrix Asub of A that is
// located col sub-matrices to the right and row sub-matrices down
// from the upper-left corner of A
template <Integer Layout, typename M>
__device__ std::enable_if_t<Layout == 1, M> GetSubMatrix(
   M A, int row, int col) {
   M Asub;
   Asub.width = BLOCK_SIZE;
   Asub.height = BLOCK_SIZE;
   Asub.stride = A.stride;
   Asub.elements = &A.elements[A.stride * BLOCK_SIZE * row + BLOCK_SIZE * col];
   return Asub;
}

//BLOCK GEMM KERNELS
//

// irregular block gemm kernel, default Layout left (L = 0)
template <typename M, typename L, typename T, typename S, typename I,
   Integer Layout = 0>
__device__ void iBlockGemmDevice(M Astart, L Bstart, L Cstart, const S rowScan,
   const S colScan, const S allscan, const T *xGrid, const T *yGrid,
   const Integer matrix_block_id, const I kWidth, const I matrixRowSize) {

   int blockid = blockIdx.x;

   // 2D block coordinates of each block id
   int blockRow = xGrid[blockid];
   int blockCol = yGrid[blockid];

   // Each thread block computes one sub-matrix Csub of C
   L Csub = GetSubMatrix<Layout>(Cstart, blockRow, blockCol);

   using CScalar = typename L::Scalar;
   // Each thread computes one element of Csub
   // by accumulating results into Cvalue
   CScalar Cvalue = 0;

   auto rs = kWidth % BLOCK_SIZE;
   auto subMatrixSize = kWidth / BLOCK_SIZE;

   if(rs > 0)
      ++subMatrixSize;

   // Thread row and column within Csub
   int row = threadIdx.x;
   int col = threadIdx.y;

   auto lastBlockRow =
      rowScan(matrix_block_id + 1) - rowScan(matrix_block_id) - 1;
   auto lastBlockCol =
      colScan(matrix_block_id + 1) - colScan(matrix_block_id) - 1;

   auto Axblock = BLOCK_SIZE;
   auto Ayblock = BLOCK_SIZE;

   /* auto col_size = Cstart.width; */
   // handle the last row block which can be less than block size.
   if(lastBlockRow == blockRow)
      Axblock = matrixRowSize - BLOCK_SIZE * lastBlockRow;

   // handle the last col block which can be less than block size.
   if(lastBlockCol == blockCol)
      Ayblock = Cstart.width - BLOCK_SIZE * lastBlockCol;

   // Loop over all the sub-matrices of A and B that are rquired to compute Csub
   // Multiply each pair of sub-matrices together and accumulate the results
   for(int m = 0; m < subMatrixSize; ++m)
   {
      auto ABxyBlock = BLOCK_SIZE;
      /* handle the last col block which can be less than block size. */
      if(m == subMatrixSize - 1 && rs > 0)
         ABxyBlock = Bstart.block_stride - BLOCK_SIZE * m;

      using AScalar = typename M::Scalar;
      using BScalar = typename L::Scalar;

      /* Shared memory used to store Asub and Bsub respectively */
      /* BLOCK_SIZE +  1 to reduce shared memory bank conflicts when column
       * major */
      __shared__ AScalar As[BLOCK_SIZE][BLOCK_SIZE + 1];
      __shared__ AScalar Bs1[BLOCK_SIZE][BLOCK_SIZE + 1];
      __shared__ AScalar Bs2[BLOCK_SIZE][BLOCK_SIZE + 1];

      /* Get sub - matrix Asub of A */
      M Asub = GetSubMatrix<Layout>(Astart, blockRow, m);

      /* Get sub - matrix Bsub of B */
      L Bsub = GetSubMatrix<Layout>(Bstart, m, blockCol);

      /* Load Asub and Bsub from device memory to shared memory
          Each thread loads one element of each sub-matrix */

      // Handle the last row block
      if(col < ABxyBlock && row < Axblock)
      { As[row][col] = GetElement<Layout>(Asub, row, col); }

      // Handle the last col block
      if(row < ABxyBlock && col < Ayblock)
      {
         auto bsub = GetElement<Layout>(Bsub, row, col);
         Bs1[row][col] = bsub.real();
         Bs2[row][col] = bsub.imag();
      }

      /* Synchronize to make sure the sub matrices are loaded
      to shared memory before starting the computation */
      __syncthreads();
      /* Multiply Asub and Bsub together */
      // handle the last col and last row block
      if(row < Axblock && col < Ayblock)
      {
         for(int e = 0; e < ABxyBlock; ++e)
         { Cvalue += As[row][e] * CScalar(Bs1[e][col], Bs2[e][col]); }
      }

      /* Synchronize to make sure that the preceding computation is done
      before loading two new sub-matrices of A and B in the next iteration */
      __syncthreads();
   }

   // handle the last col and last row block
   if(row < Axblock && col < Ayblock)
   {
      // Write Csub to device memory
      // Each thread writes one element
      SetElement<Layout>(Csub, row, col, Cvalue);
   }
}

// HOST UTILITIES

template <typename V, typename M>
void generate_block_cluster(
   V xGrid, V yGrid, M hostRowScan, M hostColScan, M hostAllScan) {

   auto slowSize = hostRowScan.extent(0) - 1;

   auto hxGrid = Kokkos::create_mirror_view(xGrid);
   auto hyGrid = Kokkos::create_mirror_view(yGrid);

   // Each block is assigned a 2D coordinate  in the matrix
   for(int i = 0; i < slowSize; i++)
   {
      auto rowBlocks = hostRowScan(i + 1) - hostRowScan(i);
      auto colBlocks = hostColScan(i + 1) - hostColScan(i);
      /* Create CUDA Dim Grid (block cluster) as it is not regular.
       * In such way to preserve coalescing by processing column wise*/
      for(int l = 0; l < colBlocks; l++)
      {
         for(int k = 0; k < rowBlocks; k++)
         {
            auto index = hostAllScan(i) + l * rowBlocks + k;
            hxGrid(index) = k;
            hyGrid(index) = l;
         }
      }
   }

   Kokkos::deep_copy(xGrid, hxGrid);
   Kokkos::deep_copy(yGrid, hyGrid);
}

// irregular block gemm kernel, default Layout left (L = 0), horizontal mops
// matrix
template <typename M, typename L, typename T, typename S, Integer Layout = 0>
__global__ void iBlockGemmKernelProjector(M A, L B, L C, const S rowScan,
   const S colScan, const S scan, const S allscan, const T *xGrid,
   const T *yGrid) {
   int blockid = blockIdx.x;

   // The id of the matrix block within the large matrix
   auto matrix_block_id = binary_search_range(allscan, blockid);
   // The start address of the A & B matrices
   auto matrix_block_row_start = scan(matrix_block_id);
   // The size of the block matrix with matrix block id.
   auto BRowSize = scan(matrix_block_id + 1) - matrix_block_row_start;
   auto CRowSize = A.height;

   M Astart = GetStartMatrix(A, matrix_block_row_start * CRowSize, CRowSize);
   /* This is a horizontal matrix. Calc. differ from the A matrix. */
   L Bstart = GetStartMatrix(B, matrix_block_row_start, BRowSize);
   // The start address of the vertical C matrix
   // matrix starts at scan(i) with size scan(i+1) - scan(i)
   L Cstart = GetStartMatrix(C, CRowSize * matrix_block_id, CRowSize, C.width);

   iBlockGemmDevice(Astart, Bstart, Cstart, rowScan, colScan, allscan, xGrid,
      yGrid, matrix_block_id, BRowSize, CRowSize);
}


// irregular block gemm kernel, default Layout left (L = 0), horizontal IN
template <typename M, typename L, typename T, typename S, Integer Layout = 0>
__global__ void iBlockGemmKernelIntegrator(M A, L B, L C, const S rowScan,
   const S colScan, const S scan, const S allscan, const T *xGrid,
   const T *yGrid) {
   int blockid = blockIdx.x;
   // The id of the matrix block within the large matrix
   auto matrix_block_id = binary_search_range(allscan, blockid);
      // The start address of the vertical C and A matrices
   auto matrix_block_row_start = scan(matrix_block_id);
   // The size of the block matrix with matrix block id.
   auto matrixRowSize = scan(matrix_block_id + 1) - matrix_block_row_start;
   // The start address of the horizontal B matrix
   auto matrix_block_col_start = matrix_block_id * B.height * C.width;

   M Astart = GetStartMatrix(A, matrix_block_row_start, matrixRowSize);
   // This is a horizontal matrix. Calc. differ from the A matrix.
   L Bstart = GetStartMatrix(B, matrix_block_col_start, B.height);
   // matrix starts at scan(i) with size scan(i+1) - scan(i)
   L Cstart = GetStartMatrix(C, matrix_block_row_start, matrixRowSize, C.width);

   iBlockGemmDevice(Astart, Bstart, Cstart, rowScan, colScan, allscan, xGrid,
      yGrid, matrix_block_id, A.width, matrixRowSize);
}

// Block Matrix multiplication - Host code - Using layout left (column major)
// matrices
template <Integer S = 0, typename M, typename L, typename V>
void iBlockGemm(const M &d_A, const L &d_B, const L &d_C, const V &rowScan,
   const V &colScan, const V &scan, const V &xGrid, const V &yGrid,
   const V &allScan, const int allTotal) {
   /* cudaEvent_t start, stop;
   cudaEventCreate(&start);
   cudaEventCreate(&stop); */

   dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
   // Invoke kernel
   /* cudaEventRecord(start); */

   if constexpr(S == 0)
   {
      iBlockGemmKernelIntegrator<<<allTotal, dimBlock>>>(d_A, d_B, d_C, rowScan,
         colScan, scan, allScan, xGrid.data(), yGrid.data());
   } else
   {
      iBlockGemmKernelProjector<<<allTotal, dimBlock>>>(d_A, d_B, d_C, rowScan,
         colScan, scan, allScan, xGrid.data(), yGrid.data());
   }

   /* cudaEventRecord(stop); */

   /* cudaEventSynchronize(stop);
   float milliseconds = 0;
   cudaEventElapsedTime(&milliseconds, start, stop); */
}

template <Integer S = 0, typename T, typename MS, typename MZ, typename V>
void applyBlockOperator(const T mspSetup, const MS &vmOps, const MZ &rOutView,
   const MZ &inView, const V &scan, const int total) {

   auto col_size = mspSetup->mult(0);
   assert(col_size == rOutView.extent(1));
   auto slowSize = mspSetup->slowSize();
   auto outRows = mspSetup->fwdSize();

   V rowScan("outRows Scan", slowSize + 1);
   V colScan("outRows Scan", slowSize + 1);
   V allScan("outRows Scan", slowSize + 1);
   auto hostRowScan = Kokkos::create_mirror_view(rowScan);
   auto hostColScan = Kokkos::create_mirror_view(colScan);
   auto hostAllScan = Kokkos::create_mirror_view(allScan);

   // build row and cols scan for each matrix using the block sizes.
   for(int i = 0; i < slowSize; i++)
   {
      if constexpr(S == 0)
      {
         auto m = mspSetup->slow(i);
         outRows = mspSetup->fast(mspSetup->fastSize(i) - 1, i) - m + 1;
      }
      auto ro = outRows % BLOCK_SIZE;
      auto rowBlocks = outRows / BLOCK_SIZE;
      auto rc = col_size % BLOCK_SIZE;
      auto colBlocks = col_size / BLOCK_SIZE;

      if(ro > 0)
         ++rowBlocks;

      if(rc > 0)
         ++colBlocks;

      hostRowScan(i + 1) = hostRowScan(i) + rowBlocks;
      hostColScan(i + 1) = hostColScan(i) + colBlocks;
      hostAllScan(i + 1) = hostAllScan(i) + rowBlocks * colBlocks;
   }

   // get totals and deep copy
   auto rowTotal = hostRowScan(slowSize);
   auto colTotal = hostColScan(slowSize);
   auto allTotal = hostAllScan(slowSize);
   Kokkos::deep_copy(rowScan, hostRowScan);
   Kokkos::deep_copy(colScan, hostColScan);
   Kokkos::deep_copy(allScan, hostAllScan);


   V xGrid("grid row index", allTotal);
   V yGrid("grid row index", allTotal);

   generate_block_cluster(xGrid, yGrid, hostRowScan, hostColScan, hostAllScan);

   // use cuda legendre data types
   using OpMatrixZC = typename CudaIALegendreOperatorTypes::OpMatrixZC;
   using OpMatrixC = typename CudaIALegendreOperatorTypes::OpMatrixC;

   // build cuda matrix types from Kokkos views
   OpMatrixC d_A;
   d_A.width = vmOps.extent(1);
   d_A.height = vmOps.extent(0);
   d_A.stride = d_A.height;
   d_A.elements = vmOps.data();

   OpMatrixZC d_B;
   d_B.width = inView.extent(1);
   d_B.height = inView.extent(0);
   d_B.stride = d_B.height;
   d_B.elements = inView.data();

   OpMatrixZC d_C;
   d_C.width = rOutView.extent(1);
   d_C.height = rOutView.extent(0);
   d_C.stride = d_C.height;
   d_C.elements = rOutView.data();

   // call the dispatch method
   iBlockGemm<S>(
      d_A, d_B, d_C, rowScan, colScan, scan, xGrid, yGrid, allScan, allTotal);
}

template <typename T, typename M>
void constantMultiplyMatrix(const T constant, const M &matrix){
   auto rows = matrix.extent(0);
   auto cols = matrix.extent(1);
   Kokkos::parallel_for(
      "constantMultiplyMatrix",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {rows, cols}),
      KOKKOS_LAMBDA(int i, int j) {
         auto result = matrix(i, j) * constant;
         matrix(i, j) = result;
      });
}

template <Integer S = 0, typename T, typename V, typename M>
void constantMultiplyMatrix(const T mspSetup, const V &scan, const M &matrix){

   using DataType = CudaIALegendreOperatorTypes::DataType;
   auto slowSize = mspSetup->slowSize();

   ViewVectorType<DataType> constants("Constant vector", slowSize);
   auto hostConstants = Kokkos::create_mirror_view(constants);

   auto sign = -1.0;
   if constexpr(S == 1){
       sign = 1.0;
   }

   for(int i = 0; i < slowSize; i++)
   {
      hostConstants(i) =
         DataType(0.0, sign * static_cast<MHDFloat>(mspSetup->slow(i)));
   }

   Kokkos::deep_copy(constants, hostConstants);

   auto rows = matrix.extent(0);
   auto cols = matrix.extent(1);

   Kokkos::parallel_for(
      "constantMultiplyMatrix using a scan",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {rows, cols}),
      KOKKOS_LAMBDA(int i, int j) {
         auto index = binary_search_range(scan, i);
         auto result = matrix(i, j) * constants(index);
         matrix(i, j) = result;
      });
}

namespace Experimental {

template <typename M, typename L, typename T, typename S, Integer Layout = 0>
__global__ void iBlockGemmExpKernel(M A, L B, L C, L C1, const S rowScan,
   const S colScan, const S scan, const S allscan, const T *xGrid,
   const T *yGrid) {

   auto col_size = C.width;

   int blockRow = blockIdx.x;
   int blockCol = blockIdx.y;

   // Thread row and column within Csub
   int row = threadIdx.x;
   int col = threadIdx.y;

   auto grow = blockDim.x * blockIdx.x + row;
   auto gcol = blockDim.y * blockIdx.y + col;

   auto matrix_block_id = binary_search_range(scan, grow);
   // Each thread block computes one sub-matrix Csub of C
   L Csub = GetSubMatrix<Layout>(C, blockRow, blockCol);

   // TODO: Check if correct. Probably yes, you use the block id from scan
   // to get the start of the in matrix. So compared to the other impl
   // this part of the B matrix can stay the same at least for now.
   // Consider doing vertical layout in.
   //
   // The start address of the horizontal B matrix
   auto matrix_block_col_start = matrix_block_id * B.height * col_size;

   using CScalar = typename L::Scalar;
   // Each thread computes one element of Csub
   // by accumulating results into Cvalue
   CScalar Cvalue = 0;
   /* cuDoubleComplex czero = make_cuDoubleComplex(0, 0);
   cuDoubleComplex CCvalue = czero; */

   auto rs = A.width % BLOCK_SIZE;
   auto subMatrixSize = A.width / BLOCK_SIZE;

   if(rs > 0)
      ++subMatrixSize;

   // Loop over all the sub-matrices of A and B that are
   // required to compute Csub
   // Multiply each pair of sub-matrices together
   // and accumulate the results
   for(int m = 0; m < subMatrixSize; ++m)
   {
      auto ABxyBlock = BLOCK_SIZE;
      /* handle the last col block which can be less than block size. */
      if(m == subMatrixSize - 1 && rs > 0)
         ABxyBlock = B.height - BLOCK_SIZE * m;

      using AScalar = typename M::Scalar;
      using BScalar = typename L::Scalar;

      /* Shared memory used to store Asub and Bsub respectively */
      /* BLOCK_SIZE +  1 to reduce shared memory bank conflicts when column
       * major */
      __shared__ AScalar As[BLOCK_SIZE][BLOCK_SIZE + 1];
      __shared__ AScalar Bs1[BLOCK_SIZE][BLOCK_SIZE + 1];
      __shared__ AScalar Bs2[BLOCK_SIZE][BLOCK_SIZE + 1];
      /* __shared__ cuDoubleComplex Bs[BLOCK_SIZE][BLOCK_SIZE + 1]; */

      /* Get sub - matrix Asub of A */
      M Asub = GetSubMatrix<Layout>(A, blockRow, m);

      // This is a horizontal matrix. Calc. differ from the A matrix.
      L Bstart = GetStartMatrix(B, matrix_block_col_start, B.height);
      /* Get sub - matrix Bsub of B */
      L Bsub = GetSubMatrix<Layout>(Bstart, m, blockCol);

      /* Load Asub and Bsub from device memory to shared memory
          Each thread loads one element of each sub-matrix */

      // Handle the last row block
      /* if(col < ABxyBlock && row < Axblock) */
      if(col < ABxyBlock && grow < C.height)
      { As[row][col] = GetElement<Layout>(Asub, row, col); }

      // Handle the last col block
      /* if(row < ABxyBlock && col < Ayblock) */
      if(row < ABxyBlock && gcol < C.width)
      {
         auto bsub = GetElement<Layout>(Bsub, row, col);
         // transpose
         /* Bs1[row_tp][col_tp] = bsub.real();
         Bs2[row_tp][col_tp] = bsub.imag(); */
         Bs1[row][col] = bsub.real();
         Bs2[row][col] = bsub.imag();
         /* Bs[row][col] = make_cuDoubleComplex(bsub.real(), bsub.imag()); */
      }

      /* Synchronize to make sure the sub matrices are loaded
      to shared memory before starting the computation */
      __syncthreads();
      /* Multiply Asub and Bsub together */

      // handle the last col and last row block
      /* if(row < Axblock && col < Ayblock) */
      if(grow < C.height && gcol < C.width)
      {
         for(int e = 0; e < ABxyBlock; ++e)
         {
            Cvalue += As[row][e] * CScalar(Bs1[e][col], Bs2[e][col]);
            /* Cvalue += As[row][e] * CScalar(Bs1[col][e], Bs2[col][e]); */
         }
      }

      /* Synchronize to make sure that the preceding computation is done
      before loading two new sub-matrices of A and B in the next iteration */
      __syncthreads();
   }

   // TODO: Get the block id of the row and block id of the col and if equal
   // store. Meaning do only comp. into your own block.
   if(grow < C.height && gcol < C.width)
   {
      /* auto val = GetElement<Layout>(Csub1, row, col); */
      SetElement<Layout>(Csub, row, col, Cvalue);

      /* auto a = GetElement<Layout>(Csub, row, col); */
      /* C.elements[grow + C.stride * gcol] = a; */
   }
}

// Block Matrix multiplication - Host code - Using layout left (column major)
// matrices calling the iblockGemmExpKernel which works on the new cuda version
// of the iblockgemm.
template <typename M, typename L, typename V>
void iBlockGemmExperimental(const M &d_A, const L &d_B, const L &d_C,
   const V &rowScan, const V &colScan, const V &scan, const V &xGrid,
   const V &yGrid, const V &allScan, const int allTotal) {
   cudaEvent_t start, stop;
   cudaEventCreate(&start);
   cudaEventCreate(&stop);

   dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

   auto ro = d_C.height % BLOCK_SIZE;
   auto rowBlocks = d_C.height / BLOCK_SIZE;
   auto rc = d_C.width % BLOCK_SIZE;
   auto colBlocks = d_C.width / BLOCK_SIZE;

   if(ro > 0)
      ++rowBlocks;

   if(rc > 0)
      ++colBlocks;

   dim3 dimGrid(rowBlocks, colBlocks);

   // Invoke kernel
   cudaEventRecord(start);

   /* iBlockCopyKernel<<<allTotal, dimBlock>>>(d_A, d_B, d_C, d_C, rowScan,
      colScan, scan, allScan, xGrid.data(), yGrid.data()); */

   iBlockGemmExpKernel<<<dimGrid, dimBlock>>>(d_A, d_B, d_C, d_C, rowScan,
      colScan, scan, allScan, xGrid.data(), yGrid.data());

   cudaEventRecord(stop);

   cudaEventSynchronize(stop);
   float milliseconds = 0;
   cudaEventElapsedTime(&milliseconds, start, stop);
}
} // namespace Experimental
#endif
} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PIALEGENDREOPERATORTYPES_HPP
