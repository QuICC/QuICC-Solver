/**
 * @file CudaIALegendreOperatorGemmUtils.hpp
 * @brief Associated Legendre based operator cuda irregular block gemm utils
 */

#ifndef QUICC_TRANSFORM_POLY_CUDAIOPERATORGEMMUTILS_HPP
#define QUICC_TRANSFORM_POLY_CUDAIOPERATORGEMMUTILS_HPP

// System includes
//

// Project includes
//

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#include "Cuda/CudaUtil.hpp"
#include "CudaIOperatorTypes.hpp"
#include "QuICC/Transform/Poly/KokkosCudaIOperatorGemmUtils.hpp"
#endif

namespace QuICC {

namespace Transform {

namespace Poly {

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
/* using DataType = cuDoubleComplex; */

// Thread block size
#define BLOCK_SIZE 16
/* #define BLOCK_SIZE 32 */

// Get a matrix element, layout left
template <Integer Layout, typename M>
__device__ std::enable_if_t<Layout == 0, typename M::Scalar> GetElement(
   const M A, int row, int col)
{
   return A.elements[col * A.stride + row];
}

// Get a matrix element, layout right
template <Integer Layout, typename M>
__device__ std::enable_if_t<Layout == 1, typename M::Scalar> GetElement(
   const M A, int row, int col)
{
   return A.elements[row * A.stride + col];
}

// Set a matrix element, layout left, cuda complex
/* template <typename M>
__device__ void SetElement(M A, int row, int col, cuDoubleComplex value) {
   A.elements[col * A.stride + row].real() = cuCreal(value);
   A.elements[col * A.stride + row].imag() = cuCimag(value);
} */

// Set a matrix element
template <Integer Layout, typename M>
__device__ std::enable_if_t<Layout == 0, void> SetElement(M A, int row, int col,
   typename M::Scalar value)
{
   A.elements[col * A.stride + row] = value;
}

// Set a matrix element
template <Integer Layout, typename M>
__device__ std::enable_if_t<Layout == 1, void> SetElement(M A, int row, int col,
   typename M::Scalar value)
{
   A.elements[row * A.stride + col] = value;
}

template <typename M>
__device__ M GetStartMatrix(M A, int offset, int stride, int width)
{
   M Asub;
   Asub.width = width;
   Asub.height = stride;
   Asub.block_stride = stride;
   Asub.stride = A.stride;
   Asub.elements = &A.elements[offset];
   return Asub;
}

template <typename M> __device__ M GetStartMatrix(M A, int offset, int stride)
{
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
__device__ std::enable_if_t<Layout == 0, M> GetSubMatrix(M A, int row, int col)
{
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
__device__ std::enable_if_t<Layout == 1, M> GetSubMatrix(M A, int row, int col)
{
   M Asub;
   Asub.width = BLOCK_SIZE;
   Asub.height = BLOCK_SIZE;
   Asub.stride = A.stride;
   Asub.elements = &A.elements[A.stride * BLOCK_SIZE * row + BLOCK_SIZE * col];
   return Asub;
}

struct Abs2Complex
{
   template <typename C> __device__ auto operator()(C v) const
   {
      return (v.real() * v.real()) + (v.imag() * v.imag());
   }
};


template <typename V, typename M, typename L, typename T, Integer Layout = 0>
__global__ void MatMulKernel(M A, L B, T C)
{
   // Block row and column
   int blockRow = blockIdx.y;
   int blockCol = blockIdx.x;

   // Each thread block computes one sub-matrix Csub of C
   T Csub = GetSubMatrix<Layout>(C, blockRow, blockCol);

   // Each thread computes one element of Csub
   // by accumulating results into Cvalue
   V Cvalue = 0;

   // Thread row and column within Csub
   int row = threadIdx.y;
   int col = threadIdx.x;

   // Loop over all the sub-matrices of A and B that are
   // required to compute Csub
   // Multiply each pair of sub-matrices together
   // and accumulate the results
   for (int m = 0; m < (A.width / BLOCK_SIZE); ++m)
   {

      // Get sub-matrix Asub of A
      M Asub = GetSubMatrix<Layout>(A, blockRow, m);

      // Get sub-matrix Bsub of B
      L Bsub = GetSubMatrix<Layout>(B, m, blockCol);

      // Shared memory used to store Asub and Bsub respectively
      __shared__ V As[BLOCK_SIZE][BLOCK_SIZE];
      __shared__ V Bs[BLOCK_SIZE][BLOCK_SIZE];

      // Load Asub and Bsub from device memory to shared memory
      // Each thread loads one element of each sub-matrix
      As[row][col] = GetElement<Layout>(Asub, row, col);
      Bs[row][col] = GetElement<Layout>(Bsub, row, col);

      // Synchronize to make sure the sub-matrices are loaded
      // before starting the computation
      __syncthreads();
      // Multiply Asub and Bsub together
      for (int e = 0; e < BLOCK_SIZE; ++e)
      {
         Cvalue += As[row][e] * Bs[e][col];
      }

      // Synchronize to make sure that the preceding
      // computation is done before loading two new
      // sub-matrices of A and B in the next iteration
      __syncthreads();
   }

   // Write Csub to device memory
   // Each thread writes one element
   SetElement<Layout>(Csub, row, col, Cvalue);
}

// Matrix multiplication - Host code
// Matrix dimensions are assumed to be multiples of BLOCK_SIZE
template <typename V, typename M, typename L, typename T, Integer Layout = 0>
void MatMul(const M A, const L B, T C)
{
   // Load A and B to device memory
   using OpMatrixC = typename CudaIOperatorTypes::OpMatrixC;

   OpMatrixC d_A;
   d_A.width = A.extent(1);
   d_A.height = A.extent(0);
   d_A.stride = d_A.height;
   d_A.elements = A.data();

   OpMatrixC d_B;
   d_B.width = B.extent(1);
   d_B.height = B.extent(0);
   d_B.stride = d_B.height;
   d_B.elements = B.data();

   OpMatrixC d_C;
   d_C.width = C.extent(1);
   d_C.height = C.extent(0);
   d_C.stride = d_C.height;
   d_C.elements = C.data();


   // Invoke kernel
   dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
   dim3 dimGrid(d_B.width / dimBlock.x, d_A.height / dimBlock.y);
   MatMulKernel<V><<<dimGrid, dimBlock>>>(d_A, d_B, d_C);
}

// BLOCK GEMM KERNELS
//

// irregular block gemm kernel, default Layout left (L = 0)
template <typename M, typename L, typename V, typename T, typename S,
   Integer Layout = 0, typename F = nullptr_t>
__device__ void iBlockGemmDevice(const F& f, M Astart, L Bstart, V Cstart,
   const S rowScan, const S colScan, const T* xGrid, const T* yGrid,
   const Integer matrix_block_id)
{

   int blockid = blockIdx.x;

   // 2D block coordinates of each block id
   int blockRow = xGrid[blockid];
   int blockCol = yGrid[blockid];

   // Each thread block computes one sub-matrix Csub of C
   V Csub = GetSubMatrix<Layout>(Cstart, blockRow, blockCol);

   using RScalar = typename V::Scalar;
   using CScalar = typename L::Scalar;
   // Each thread computes one element of Csub
   // by accumulating results into Cvalue
   CScalar Cvalue = 0;

   auto rs = Bstart.height % BLOCK_SIZE;
   auto subMatrixSize = Bstart.height / BLOCK_SIZE;

   if (rs > 0)
   {
      ++subMatrixSize;
   }

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
   if (lastBlockRow == blockRow)
   {
      Axblock = Cstart.height - BLOCK_SIZE * lastBlockRow;
   }

   // handle the last col block which can be less than block size.
   if (lastBlockCol == blockCol)
   {
      Ayblock = Cstart.width - BLOCK_SIZE * lastBlockCol;
   }

   // Loop over all the sub-matrices of A and B that are rquired to compute Csub
   // Multiply each pair of sub-matrices together and accumulate the results
   for (int m = 0; m < subMatrixSize; ++m)
   {
      auto ABxyBlock = BLOCK_SIZE;
      /* handle the last col block which can be less than block size. */
      if (m == subMatrixSize - 1 && rs > 0)
      {
         ABxyBlock = Bstart.block_stride - BLOCK_SIZE * m;
      }

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
      if (col < ABxyBlock && row < Axblock)
      {
         As[row][col] = GetElement<Layout>(Asub, row, col);
      }

      // Handle the last col block
      if (row < ABxyBlock && col < Ayblock)
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
      if (row < Axblock && col < Ayblock)
      {
         for (int e = 0; e < ABxyBlock; ++e)
         {
            Cvalue += As[row][e] * CScalar(Bs1[e][col], Bs2[e][col]);
         }
      }

      /* Synchronize to make sure that the preceding computation is done
      before loading two new sub-matrices of A and B in the next iteration */
      __syncthreads();
   }

   // handle the last col and last row block
   if (row < Axblock && col < Ayblock)
   {
      RScalar result = 0;
      if constexpr (!std::is_same_v<F, nullptr_t>)
      {
         result = f(Cvalue);
      }
      else
      {
         static_assert(std::is_same_v<V, L>);
         result = Cvalue;
      }
      // Write Csub to device memory
      // Each thread writes one element
      SetElement<Layout>(Csub, row, col, result);
   }
}

// HOST UTILITIES
//
//


// irregular block gemm kernel, default Layout left (L = 0), horizontal IN
template <typename M, typename L, typename V, typename T, typename S,
   Integer Layout = 0, typename F = nullptr_t>
__global__ void iBlockGemmKernelWorlandProjector(M A, L B, V C, const S rowScan,
   const S colScan, const S scan, const S allscan, const T* xGrid,
   const T* yGrid, const F& f)
{
   int blockid = blockIdx.x;
   // The id of the matrix block within the large matrix
   auto matrix_block_id = binary_search_range(allscan, blockid);
   // The start address of the vertical C and A matrices
   auto matrix_block_col_start = scan(matrix_block_id);
   // The size of the block matrix with matrix block id.
   auto matrixColSize = scan(matrix_block_id + 1) - matrix_block_col_start;
   // The start address of the horizontal B matrix
   auto matrix_block_row_start = matrix_block_id * C.height * B.height;

   M Astart = GetStartMatrix(A, matrix_block_row_start, C.height);
   // This is a horizontal matrix. Calc. differ from the A matrix.
   L Bstart = GetStartMatrix(B, matrix_block_col_start * B.height, B.height);
   // matrix starts at scan(i) with size scan(i+1) - scan(i)
   V Cstart = GetStartMatrix(C, matrix_block_col_start * C.height, C.height,
      matrixColSize);

   iBlockGemmDevice(f, Astart, Bstart, Cstart, rowScan, colScan, xGrid, yGrid,
      matrix_block_id);
}

// irregular block gemm kernel, default Layout left (L = 0), horizontal mops
// matrix
template <typename M, typename L, typename V, typename T, typename S,
   Integer Layout = 0, typename F = nullptr_t>
__global__ void iBlockGemmKernelProjector(M A, L B, V C, const S rowScan,
   const S colScan, const S scan, const S allscan, const T* xGrid,
   const T* yGrid, const F& f)
{
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
   V Cstart = GetStartMatrix(C, CRowSize * matrix_block_id, CRowSize, C.width);

   iBlockGemmDevice(f, Astart, Bstart, Cstart, rowScan, colScan, xGrid, yGrid,
      matrix_block_id);
}

// irregular block gemm kernel, default Layout left (L = 0), horizontal IN
template <typename M, typename L, typename V, typename T, typename S,
   Integer Layout = 0, typename F = nullptr_t>
__global__ void iBlockGemmKernelReductor(M A, L B, V C, const S rowScan,
   const S colScan, const S scan, const S allscan, const T* xGrid,
   const T* yGrid, const F& f)
{
   int blockid = blockIdx.x;
   // The id of the matrix block within the large matrix
   auto matrix_block_id = binary_search_range(allscan, blockid);
   // The start address of the vertical C and A matrices
   auto matrix_block_row_start = scan(matrix_block_id);
   // The size of the block matrix with matrix block id.
   auto matrixRowSize = scan(matrix_block_id + 1) - matrix_block_row_start;
   // The start address of the horizontal B matrix
   auto matrix_block_col_start = matrix_block_id * A.width;

   M Astart = GetStartMatrix(A, matrix_block_row_start, matrixRowSize);
   // This is a horizontal matrix. Calc. differ from the A matrix.
   L Bstart = GetStartMatrix(B, matrix_block_col_start, A.width);
   // matrix starts at scan(i) with size scan(i+1) - scan(i)
   V Cstart = GetStartMatrix(C, matrix_block_row_start, matrixRowSize, C.width);

   iBlockGemmDevice(f, Astart, Bstart, Cstart, rowScan, colScan, xGrid, yGrid,
      matrix_block_id);
}


// irregular block gemm kernel, default Layout left (L = 0), horizontal IN
template <typename M, typename L, typename V, typename T, typename S,
   Integer Layout = 0, typename F = nullptr_t>
__global__ void iBlockGemmKernelWorlandIntegrator(M A, L B, V C,
   const S rowScan, const S colScan, const S scan, const S allscan,
   const T* xGrid, const T* yGrid, const F& f)
{
   int blockid = blockIdx.x;
   // The id of the matrix block within the large matrix
   auto matrix_block_id = binary_search_range(allscan, blockid);
   // The start address of the vertical C and A matrices
   auto matrix_block_col_start = scan(matrix_block_id);
   // The size of the block matrix with matrix block id.
   auto matrixColSize = scan(matrix_block_id + 1) - matrix_block_col_start;
   // The start address of the horizontal B matrix
   auto matrix_block_row_start = matrix_block_id * C.height;

   M Astart = GetStartMatrix(A, matrix_block_row_start, C.height);
   // This is a horizontal matrix. Calc. differ from the A matrix.
   L Bstart = GetStartMatrix(B, matrix_block_col_start * B.height, B.height);
   // matrix starts at scan(i) with size scan(i+1) - scan(i)
   V Cstart = GetStartMatrix(C, matrix_block_col_start * C.height, C.height,
      matrixColSize);

   iBlockGemmDevice(f, Astart, Bstart, Cstart, rowScan, colScan, xGrid, yGrid,
      matrix_block_id);
}


// irregular block gemm kernel, default Layout left (L = 0), horizontal IN
template <typename M, typename L, typename V, typename T, typename S,
   Integer Layout = 0, typename F = nullptr_t>
__global__ void iBlockGemmKernelIntegrator(M A, L B, V C, const S rowScan,
   const S colScan, const S scan, const S allscan, const T* xGrid,
   const T* yGrid, const F& f)
{
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
   V Cstart = GetStartMatrix(C, matrix_block_row_start, matrixRowSize, C.width);

   iBlockGemmDevice(f, Astart, Bstart, Cstart, rowScan, colScan, xGrid, yGrid,
      matrix_block_id);
}

// Block Matrix multiplication - Host code - Using layout left (column major)
// matrices
template <Integer S = 0, typename M, typename L, typename Z, typename V,
   typename F = nullptr_t>
void iBlockGemm(const M& d_A, const L& d_B, const Z& d_C, const V& rowScan,
   const V& colScan, const V& scan, const V& xGrid, const V& yGrid,
   const V& allScan, const int allTotal, const F& f)
{

   dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

   if constexpr (S == 0)
   {
      iBlockGemmKernelIntegrator<<<allTotal, dimBlock>>>(d_A, d_B, d_C, rowScan,
         colScan, scan, allScan, xGrid.data(), yGrid.data(), f);
   }
   else if (S == 1)
   {
      iBlockGemmKernelProjector<<<allTotal, dimBlock>>>(d_A, d_B, d_C, rowScan,
         colScan, scan, allScan, xGrid.data(), yGrid.data(), f);
   }
   else if (S == 2)
   {
      iBlockGemmKernelReductor<<<allTotal, dimBlock>>>(d_A, d_B, d_C, rowScan,
         colScan, scan, allScan, xGrid.data(), yGrid.data(), f);
   }
   else if (S == 3)
   {
      iBlockGemmKernelWorlandIntegrator<<<allTotal, dimBlock>>>(d_A, d_B, d_C,
         rowScan, colScan, scan, allScan, xGrid.data(), yGrid.data(), f);
   }
   else if (S == 4)
   {
      iBlockGemmKernelWorlandProjector<<<allTotal, dimBlock>>>(d_A, d_B, d_C,
         rowScan, colScan, scan, allScan, xGrid.data(), yGrid.data(), f);
   }


   cudaErrChk(cudaGetLastError());
   cudaErrChk(cudaDeviceSynchronize());
}

template <Integer S = 0, typename T, typename MS, typename MV, typename MZ,
   typename V, typename F = nullptr_t>
void applyBlockOperator(const T mspSetup, const MS& vmOps, const MV& rOutView,
   const MZ& inView, const V& scan, const int total, const F& f = F())
{

   auto slowSize = mspSetup->slowSize();
   auto outRows = mspSetup->fwdSize();

   V rowScan("outRows Scan", slowSize + 1);
   V colScan("outRows Scan", slowSize + 1);
   V allScan("outRows Scan", slowSize + 1);
   auto hostRowScan = Kokkos::create_mirror_view(rowScan);
   auto hostColScan = Kokkos::create_mirror_view(colScan);
   auto hostAllScan = Kokkos::create_mirror_view(allScan);

   // build row and cols scan for each matrix using the block sizes.
   for (int i = 0; i < slowSize; i++)
   {
      if constexpr (S != 1 && S != 4)
      {
         outRows = mspSetup->fastSize(i);
      }
      auto col_size = mspSetup->mult(i);

      auto ro = outRows % BLOCK_SIZE;
      auto rowBlocks = outRows / BLOCK_SIZE;
      auto rc = col_size % BLOCK_SIZE;
      auto colBlocks = col_size / BLOCK_SIZE;

      if (ro > 0)
      {
         ++rowBlocks;
      }

      if (rc > 0)
      {
         ++colBlocks;
      }

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

   // use cuda poly data types
   using OpMatrixA = typename CudaIOperatorTypes::OpMatrixC;
   using OpMatrixB = typename CudaIOperatorTypes::OpMatrixZC;


   // build cuda matrix types from Kokkos views
   OpMatrixA d_A;
   d_A.width = vmOps.extent(1);
   d_A.height = vmOps.extent(0);
   d_A.stride = d_A.height;
   d_A.elements = vmOps.data();

   OpMatrixB d_B;
   d_B.width = inView.extent(1);
   d_B.height = inView.extent(0);
   d_B.stride = d_B.height;
   d_B.elements = inView.data();

   if constexpr (!std::is_same_v<F, nullptr_t>)
   {
      using OpMatrixCC = typename CudaIOperatorTypes::OpMatrixC;

      OpMatrixCC d_C;
      d_C.width = rOutView.extent(1);
      d_C.height = rOutView.extent(0);
      d_C.stride = d_C.height;
      d_C.elements = rOutView.data();
      // call the
      iBlockGemm<S>(d_A, d_B, d_C, rowScan, colScan, scan, xGrid, yGrid,
         allScan, allTotal, f);
   }
   else
   {

      using OpMatrixCC = typename CudaIOperatorTypes::OpMatrixZC;
      OpMatrixCC d_C;
      d_C.width = rOutView.extent(1);
      d_C.height = rOutView.extent(0);
      d_C.stride = d_C.height;
      d_C.elements = rOutView.data();
      // call the
      iBlockGemm<S>(d_A, d_B, d_C, rowScan, colScan, scan, xGrid, yGrid,
         allScan, allTotal, f);
   }
}

template <typename T, typename M>
void constantMultiplyMatrix(const T constant, const M& matrix)
{

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
void constantMultiplyMatrix(const T mspSetup, const V& scan, const M& matrix)
{

   using DataType = CudaIOperatorTypes::DataType;
   auto slowSize = mspSetup->slowSize();

   ViewVectorType<DataType> constants("Constant vector", slowSize);
   auto hostConstants = Kokkos::create_mirror_view(constants);

   auto sign = -1.0;
   if constexpr (S == 1)
   {
      sign = 1.0;
   }

   for (int i = 0; i < slowSize; i++)
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
__global__ void iBlockGemmExpKernel(M A, L B, L C, const S rowScan,
   const S colScan, const S scan, const S allscan, const T* xGrid,
   const T* yGrid)
{

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

   auto rs = A.width % BLOCK_SIZE;
   auto subMatrixSize = A.width / BLOCK_SIZE;

   if (rs > 0)
   {
      ++subMatrixSize;
   }

   // Loop over all the sub-matrices of A and B that are
   // required to compute Csub
   // Multiply each pair of sub-matrices together
   // and accumulate the results
   for (int m = 0; m < subMatrixSize; ++m)
   {
      auto ABxyBlock = BLOCK_SIZE;
      /* handle the last col block which can be less than block size. */
      if (m == subMatrixSize - 1 && rs > 0)
      {
         ABxyBlock = B.height - BLOCK_SIZE * m;
      }

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
      if (col < ABxyBlock && grow < C.height)
      {
         As[row][col] = GetElement<Layout>(Asub, row, col);
      }

      // Handle the last col block
      /* if(row < ABxyBlock && col < Ayblock) */
      if (row < ABxyBlock && gcol < C.width)
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
      if (grow < C.height && gcol < C.width)
      {
         for (int e = 0; e < ABxyBlock; ++e)
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
   if (grow < C.height && gcol < C.width)
   {
      SetElement<Layout>(Csub, row, col, Cvalue);
   }
}

// Block Matrix multiplication - Host code - Using layout left (column major)
// matrices calling the iblockGemmExpKernel which works on the new cuda version
// of the iblockgemm.
template <typename M, typename L, typename V>
void iBlockGemmExperimental(const M& d_A, const L& d_B, const L& d_C,
   const V& rowScan, const V& colScan, const V& scan, const V& xGrid,
   const V& yGrid, const V& allScan, const int allTotal)
{

   dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

   auto ro = d_C.height % BLOCK_SIZE;
   auto rowBlocks = d_C.height / BLOCK_SIZE;
   auto rc = d_C.width % BLOCK_SIZE;
   auto colBlocks = d_C.width / BLOCK_SIZE;

   if (ro > 0)
   {
      ++rowBlocks;
   }

   if (rc > 0)
   {
      ++colBlocks;
   }

   dim3 dimGrid(rowBlocks, colBlocks);

   // Invoke kernel
   iBlockGemmExpKernel<<<dimGrid, dimBlock>>>(d_A, d_B, d_C, rowScan, colScan,
      scan, allScan, xGrid.data(), yGrid.data());
}
} // namespace Experimental
#endif
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_PIOPERATORTYPES_HPP
