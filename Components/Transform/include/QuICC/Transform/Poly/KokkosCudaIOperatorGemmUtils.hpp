/**
 * @file PIOperatorUtils.hpp
 * @brief Poly based operator kokkos irregular block gemm utils
 */

#ifndef QUICC_TRANSFORM_POLY_PIOPERATORGEMMUTILS_HPP
#define QUICC_TRANSFORM_POLY_PIOPERATORGEMMUTILS_HPP

// System includes
//

// Project includes
//

#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/KokkosUtils.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

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


// Requires rOut to be a vertical matrix instead of horizontal.
// Better coalescing achieved with very good occupancy.
// However requires specialized copy of rout into its original format
template <typename R, typename L, typename T, typename V>
struct ApplyUnitOperator {

   ApplyUnitOperator(R rv, L iv, T vo, V s, int ir, int cs)
       : rOutView(rv), inView(iv), Ops(vo), scan(s), inRows(ir), col_size(cs) {}

   KOKKOS_INLINE_FUNCTION
   void operator()(const int row, const int col) const {
      auto index = binary_search_range(scan, row);
      auto colin = col + index * col_size;

      for(int j = 0; j < inRows; j++)
      { rOutView(row, col) += Ops(row, j) * inView(j, colin); }
   }

   R rOutView;
   L inView;
   T Ops;
   V scan;
   int inRows;
   int col_size;
};

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

   using DataType =  typename KokkosIOperatorTypes::DataType;
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

/* using DataType = cuDoubleComplex; */

// GEMM UTILS
//
// Thread block size
#define BLOCK_SIZE 16
/* #define BLOCK_SIZE 32 */


// Choose Iteration Layout for copying data from global memory into scratch
// On CPUs it is more important to have consecutive write,
// On GPUs it is more important to not jump around in global memory, i.e. have
// coallesced loads
template <class ExecSpace, class LayoutA, class LayoutAScratch>
struct impl_gemm_choose_copy_layout {
   using type = LayoutAScratch;
};

#ifdef KOKKOS_ENABLE_CUDA
template <class LayoutA, class LayoutAScratch>
struct impl_gemm_choose_copy_layout<Kokkos::Cuda, LayoutA, LayoutAScratch> {
   using type = LayoutA;
};
#endif

#ifdef KOKKOS_ENABLE_HIP
template <class LayoutA, class LayoutAScratch>
struct impl_gemm_choose_copy_layout<Kokkos::Experimental::HIP, LayoutA,
   LayoutAScratch> {
   using type = LayoutA;
};
#endif

// DeepCopy matrix block into scratch
template <class TeamHandle, class ViewTypeScratch, class ViewType, class Layout,
   int blockDim_i, int blockDim_j, int Transpose>
struct impl_deep_copy_matrix_block;

template <class TeamHandle, class ViewTypeScratch, class ViewType, class Layout,
   int blockDim_i, int blockDim_j>
struct impl_deep_copy_matrix_block<TeamHandle, ViewTypeScratch, ViewType,
   Layout, blockDim_i, blockDim_j, 0> {
   typedef typename ViewType::non_const_value_type value_type;
   typedef Kokkos::Details::ArithTraits<value_type> ATV;

   KOKKOS_INLINE_FUNCTION
   static void copy(const TeamHandle &team, const ViewTypeScratch &A_scr,
      const ViewType &A, const int &offset_i, const int &offset_j) {
      if(offset_i + blockDim_i <= A.extent_int(0) &&
         offset_j + blockDim_j <= A.extent_int(1))
      {
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, blockDim_j), [&](const int j) {
#ifndef KOKKOS_IMPL_BATCHED_GEMM_GCC_CXX14_WORKAROUND
               const int idx_j = offset_j + j;
#endif
               Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, blockDim_i),
                  [&](const int i) {
#ifdef KOKKOS_IMPL_BATCHED_GEMM_GCC_CXX14_WORKAROUND
                     const int idx_j = offset_j + j;
#endif
                     const int idx_i = offset_i + i;
                     A_scr(i, j) = A(idx_i, idx_j);
                  });
            });
      } else
      {
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, blockDim_j), [&](const int j) {
#ifndef KOKKOS_IMPL_BATCHED_GEMM_GCC_CXX14_WORKAROUND
               int idx_j = offset_j + j;
#endif
               Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, blockDim_i),
                  [&](const int i) {
#ifdef KOKKOS_IMPL_BATCHED_GEMM_GCC_CXX14_WORKAROUND
                     int idx_j = offset_j + j;
#endif
                     const int idx_i = offset_i + i;
                     A_scr(i, j) =
                        idx_i < A.extent_int(0) && idx_j < A.extent_int(1)
                           ? A(idx_i, idx_j)
                           : ATV::zero();
                  });
            });
      }
   }
};

template <class TeamHandle, class ViewTypeScratch, class ViewType,
   int blockDim_i, int blockDim_j>
struct impl_deep_copy_matrix_block<TeamHandle, ViewTypeScratch, ViewType,
   Kokkos::LayoutRight, blockDim_i, blockDim_j, 0> {
   typedef typename ViewType::non_const_value_type value_type;
   typedef Kokkos::Details::ArithTraits<value_type> ATV;

   KOKKOS_INLINE_FUNCTION
   static void copy(const TeamHandle &team, const ViewTypeScratch &A_scr_real, const ViewTypeScratch &A_scr_imag,
      const ViewType &A, const int &offset_i, const int &offset_j) {
      if(offset_i + blockDim_i <= A.extent_int(0) &&
         offset_j + blockDim_j <= A.extent_int(1))
      {
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, blockDim_i), [&](const int i) {
               const int idx_i = offset_i + i;
               Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, blockDim_j),
                  [&](const int j) {
                     const int idx_j = offset_j + j;
                     auto complex_value = A(idx_i, idx_j);
                     A_scr_real(i, j) = complex_value.real();
                     A_scr_imag(i, j) = complex_value.imag();
                  });
            });
      } else
      {
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, blockDim_i), [&](const int i) {
#ifndef KOKKOS_IMPL_BATCHED_GEMM_GCC_CXX14_WORKAROUND
               int idx_i = offset_i + i;
#endif
               Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, blockDim_j),
                  [&](const int j) {
#ifdef KOKKOS_IMPL_BATCHED_GEMM_GCC_CXX14_WORKAROUND
                     int idx_i = offset_i + i;
#endif
                     const int idx_j = offset_j + j;
                     auto complex_value =
                        idx_i < A.extent_int(0) && idx_j < A.extent_int(1)
                           ? A(idx_i, idx_j)
                           : ATV::zero();
                     A_scr_real(i, j) = complex_value.real();
                     A_scr_imag(i, j) = complex_value.imag();
                  });
            });
      }
   }
};


template <class TeamHandle, class ViewType, class ViewTypeScratch, class Layout,
   int blockDim_i, int blockDim_j>
struct impl_update_matrix_block {
   typedef typename ViewType::non_const_value_type value_type;
   typedef Kokkos::Details::ArithTraits<value_type> ATV;

   KOKKOS_INLINE_FUNCTION
   static void update(const TeamHandle &team, const value_type &beta,
      const ViewType &A, const value_type &alpha, const ViewTypeScratch &A_scr,
      const int &offset_i, const int &offset_j) {
      if(offset_i + blockDim_i <= A.extent_int(0) &&
         offset_j + blockDim_j <= A.extent_int(1))
      {
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, blockDim_j), [&](const int j) {
               const int idx_j = offset_j + j;
               if(beta == ATV::zero())
               {
                  Kokkos::parallel_for(
                     Kokkos::ThreadVectorRange(team, blockDim_i),
                     [&](const int i) {
                        const int idx_i = offset_i + i;
                        A(idx_i, idx_j) = alpha * A_scr(i, j);
                     });
               } else
               {
                  Kokkos::parallel_for(
                     Kokkos::ThreadVectorRange(team, blockDim_i),
                     [&](const int i) {
                        const int idx_i = offset_i + i;
                        A(idx_i, idx_j) =
                           beta * A(idx_i, idx_j) + alpha * A_scr(i, j);
                     });
               }
            });
      } else
      {
         const int range_i = offset_i + blockDim_i <= A.extent_int(0)
                                ? blockDim_i
                                : A.extent_int(0) % blockDim_i;
         const int range_j = offset_j + blockDim_j <= A.extent_int(1)
                                ? blockDim_j
                                : A.extent_int(1) % blockDim_j;
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, range_j), [&](const int j) {
               const int idx_j = offset_j + j;
               if(beta == ATV::zero())
               {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, range_i),
                     [&](const int i) {
                        const int idx_i = offset_i + i;
                        A(idx_i, idx_j) = alpha * A_scr(i, j);
                     });
               } else
               {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, range_i),
                     [&](const int i) {
                        const int idx_i = offset_i + i;
                        A(idx_i, idx_j) =
                           beta * A(idx_i, idx_j) + alpha * A_scr(i, j);
                     });
               }
            });
      }
   }
};

template <class TeamHandle, class ViewType, class ViewTypeScratch,
   int blockDim_i, int blockDim_j>
struct impl_update_matrix_block<TeamHandle, ViewType, ViewTypeScratch,
   Kokkos::LayoutRight, blockDim_i, blockDim_j> {
   typedef typename ViewType::non_const_value_type value_type;
   typedef Kokkos::Details::ArithTraits<value_type> ATV;

   KOKKOS_INLINE_FUNCTION
   static void update(const TeamHandle &team, const value_type &beta,
      const ViewType &A, const value_type &alpha, const ViewTypeScratch &A_scr,
      const int &offset_i, const int &offset_j) {
      if(offset_i + blockDim_i <= A.extent_int(0) &&
         offset_j + blockDim_j <= A.extent_int(1))
      {
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, blockDim_i), [&](const int i) {
               const int idx_i = offset_i + i;
               if(beta == ATV::zero())
               {
                  Kokkos::parallel_for(
                     Kokkos::ThreadVectorRange(team, blockDim_j),
                     [&](const int j) {
                        const int idx_j = offset_j + j;
                        A(idx_i, idx_j) = alpha * A_scr(i, j);
                     });
               } else
               {
                  Kokkos::parallel_for(
                     Kokkos::ThreadVectorRange(team, blockDim_j),
                     [&](const int j) {
                        const int idx_j = offset_j + j;
                        A(idx_i, idx_j) =
                           beta * A(idx_i, idx_j) + alpha * A_scr(i, j);
                     });
               }
            });
      } else
      {
         const int range_i = offset_i + blockDim_i <= A.extent_int(0)
                                ? blockDim_i
                                : A.extent_int(0) % blockDim_i;
         const int range_j = offset_j + blockDim_j <= A.extent_int(1)
                                ? blockDim_j
                                : A.extent_int(1) % blockDim_j;
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, range_i), [&](const int i) {
               const int idx_i = offset_i + i;
               if(beta == ATV::zero())
               {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, range_j),
                     [&](const int j) {
                        const int idx_j = offset_j + j;
                        A(idx_i, idx_j) = alpha * A_scr(i, j);
                     });
               } else
               {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, range_j),
                     [&](const int j) {
                        const int idx_j = offset_j + j;
                        A(idx_i, idx_j) =
                           beta * A(idx_i, idx_j) + alpha * A_scr(i, j);
                     });
               }
            });
      }
   }
};


// Compute a single A block 8 B block, also do an in-place no-additional
// blocking team GEMM
template <class TeamHandle, class ViewTypeA, class ViewTypeB, class ViewTypeC>
KOKKOS_INLINE_FUNCTION void impl_team_gemm_block_complex(const TeamHandle &team,
   const ViewTypeC &CR, const ViewTypeC &CI, const ViewTypeA &A,
   const ViewTypeB &BR, const ViewTypeB &BI) {
   typedef typename ViewTypeA::non_const_value_type Scalar;
   typedef typename ViewTypeC::non_const_value_type ScalarC;
// GNU COMPILER BUG WORKAROUND
#if defined(KOKKOS_COMPILER_GNU) &&                                            \
   (!defined(__CUDA_ARCH__) || !defined(__HIP_DEVICE_COMPILE__))
   int blockA0 = A.extent_int(0);
   int blockA1 = A.extent_int(1);
   int blockB1 = BR.extent_int(1);
#else
   const int blockA0 = A.extent_int(0);
   const int blockA1 = A.extent_int(1);
   const int blockB1 = BR.extent_int(1);
#endif
   Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, blockA0), [&](const int i) {
#ifndef KOKKOSKERNELS_ENABLE_OMP_SIMD
         Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(team, blockB1 / 4), [&](const int B_j) {
#else
#pragma omp simd
    for(int B_j=0; B_j<blockB1/4; B_j++) {
#endif
               ScalarC C_ij0 = 0;
               ScalarC C_ij1 = 0;
               ScalarC C_ij2 = 0;
               ScalarC C_ij3 = 0;
               for(int j = 0; j < blockA1; j++)
               {
                  Scalar A_ij = A(i, j);
                  C_ij0 += A_ij * BR(j, B_j);
                  C_ij1 += A_ij * BR(j, B_j + blockB1 / 4);
                  C_ij2 += A_ij * BR(j, B_j + 2 * blockB1 / 4);
                  C_ij3 += A_ij * BR(j, B_j + 3 * blockB1 / 4);
               }

               CR(i, B_j) += C_ij0;
               CR(i, B_j + blockB1 / 4) += C_ij1;
               CR(i, B_j + 2 * blockB1 / 4) += C_ij2;
               CR(i, B_j + 3 * blockB1 / 4) += C_ij3;

               ScalarC CI_ij0 = 0;
               ScalarC CI_ij1 = 0;
               ScalarC CI_ij2 = 0;
               ScalarC CI_ij3 = 0;
               for(int j = 0; j < blockA1; j++)
               {
                  Scalar A_ij = A(i, j);
                  CI_ij0 += A_ij * BI(j, B_j);
                  CI_ij1 += A_ij * BI(j, B_j + blockB1 / 4);
                  CI_ij2 += A_ij * BI(j, B_j + 2 * blockB1 / 4);
                  CI_ij3 += A_ij * BI(j, B_j + 3 * blockB1 / 4);
               }
               CI(i, B_j) += CI_ij0;
               CI(i, B_j + blockB1 / 4) += CI_ij1;
               CI(i, B_j + 2 * blockB1 / 4) += CI_ij2;
               CI(i, B_j + 3 * blockB1 / 4) += CI_ij3;
#ifndef KOKKOSKERNELS_ENABLE_OMP_SIMD
            });
#else
    }
#endif
      });
}
// Compute a single A block 8 B block, also do an in-place no-additional
// blocking team GEMM
template <class TeamHandle, class ViewTypeA, class ViewTypeB, class ViewTypeC>
KOKKOS_INLINE_FUNCTION void impl_team_gemm_block(const TeamHandle &team,
   const ViewTypeC &C, const ViewTypeA &A, const ViewTypeB &B) {
   typedef typename ViewTypeC::non_const_value_type ScalarC;
// GNU COMPILER BUG WORKAROUND
#if defined(KOKKOS_COMPILER_GNU) &&                                            \
   (!defined(__CUDA_ARCH__) || !defined(__HIP_DEVICE_COMPILE__))
   int blockA0 = A.extent_int(0);
   int blockA1 = A.extent_int(1);
   int blockB1 = B.extent_int(1);
#else
   const int blockA0 = A.extent_int(0);
   const int blockA1 = A.extent_int(1);
   const int blockB1 = B.extent_int(1);
#endif
   Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, blockA0), [&](const int i) {
#ifndef KOKKOSKERNELS_ENABLE_OMP_SIMD
         Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(team, blockB1 / 4), [&](const int B_j) {
#else
#pragma omp simd
    for(int B_j=0; B_j<blockB1/4; B_j++) {
#endif
               ScalarC C_ij0 = 0;
               ScalarC C_ij1 = 0;
               ScalarC C_ij2 = 0;
               ScalarC C_ij3 = 0;
               for(int j = 0; j < blockA1; j++)
               {
                  ScalarC A_ij = A(i, j);
                  C_ij0 += A_ij * B(j, B_j);
                  C_ij1 += A_ij * B(j, B_j + blockB1 / 4);
                  C_ij2 += A_ij * B(j, B_j + 2 * blockB1 / 4);
                  C_ij3 += A_ij * B(j, B_j + 3 * blockB1 / 4);
               }
               C(i, B_j) += C_ij0;
               C(i, B_j + blockB1 / 4) += C_ij1;
               C(i, B_j + 2 * blockB1 / 4) += C_ij2;
               C(i, B_j + 3 * blockB1 / 4) += C_ij3;
#ifndef KOKKOSKERNELS_ENABLE_OMP_SIMD
            });
#else
    }
#endif
      });
}

template <int TransposeA, int TransposeB> struct impl_gemm_label;

template <> struct impl_gemm_label<0, 0> {
   static constexpr const char *label = "KokkosBlas::gemm[NN]";
};

template <class ExecSpace, class ViewTypeA, class ViewTypeB, class ViewTypeC,
   int blockA0, int blockA1, int blockB1, int TransposeA, int TransposeB>
struct GEMMImpl {
   ViewTypeA A;
   ViewTypeB B;
   ViewTypeC C;
   typedef typename ViewTypeA::non_const_value_type ScalarA;
   typedef typename ViewTypeB::non_const_value_type ScalarB;
   typedef typename ViewTypeC::non_const_value_type ScalarC;

   const int num_blocks_0;
   const int num_blocks_1;
   int scratch_level;

   ScalarC alpha, beta;
   typedef Kokkos::View<ScalarA[blockA0][blockA1], Kokkos::LayoutLeft,
      typename ExecSpace::scratch_memory_space>
      ViewTypeAScratch;
   typedef Kokkos::View<ScalarA[blockA1][blockB1], Kokkos::LayoutRight,
      typename ExecSpace::scratch_memory_space>
      ViewTypeBScratch;
   typedef Kokkos::View<ScalarC[blockA0][blockB1], Kokkos::LayoutRight,
      typename ExecSpace::scratch_memory_space>
      ViewTypeCScratch;

   GEMMImpl(const ViewTypeA &A_, const ViewTypeB &B_, const ViewTypeC &C_)
       : A(A_), B(B_), C(C_),
         num_blocks_0((C.extent_int(0) + blockA0 - 1) / blockA0),
         num_blocks_1((C.extent_int(1) + blockB1 - 1) / blockB1) {
      scratch_level = 0;
      alpha = 1;
      beta = 1;
   }

   GEMMImpl(const ScalarA &alpha_, const ViewTypeA &A_, const ViewTypeB &B_,
      const ScalarC &beta_, const ViewTypeC &C_)
       : A(A_), B(B_), C(C_),
         num_blocks_0((C.extent_int(0) + blockA0 - 1) / blockA0),
         num_blocks_1((C.extent_int(1) + blockB1 - 1) / blockB1) {
      scratch_level = 0;
      alpha = alpha_;
      beta = beta_;
   }

   void run(
      const ExecSpace &space, int team_size, int vector_length, int scr_level) {
      scratch_level = scr_level;
      int scratch_memory_size = ViewTypeAScratch::shmem_size() +
                                ViewTypeBScratch::shmem_size() +
                                ViewTypeCScratch::shmem_size();

#if defined(KOKKOS_ENABLE_HIP)
      // Note lbv, 10/29/20: The LaunchBounds<384,2> leads
      // to an error with HIP as the heuristics on that platform
      // yield an optimal_num_blocks=0 which means no ressources
      // are allocated... Switching to LaunchBounds<384,2> fixes
      // that problem but I'm not sure if that it a good perf
      // parameter or why it is set to 2 for Cuda?
      Kokkos::TeamPolicy<ExecSpace, Kokkos::LaunchBounds<384, 0>> policy(
         space, num_blocks_0 * num_blocks_1, team_size, vector_length);
#else
      Kokkos::TeamPolicy<ExecSpace, Kokkos::LaunchBounds<256, 2>> policy(
         space, num_blocks_0 * num_blocks_1, team_size, vector_length);
#endif

      Kokkos::parallel_for(impl_gemm_label<TransposeA, TransposeB>::label,
         policy.set_scratch_size(
            scratch_level, Kokkos::PerTeam(scratch_memory_size)),
         *this);
   }

   KOKKOS_INLINE_FUNCTION
   void operator()(
      const typename Kokkos::TeamPolicy<ExecSpace>::member_type &team) const {
      // This team is responsible for computing a single block of C
      const int blockId = team.league_rank();
      const int num_blocks = num_blocks_1;
      const int i_offset = (blockId / num_blocks) * blockA0;
      const int j_offset = (blockId % num_blocks) * blockB1;

      ViewTypeAScratch A_scr(team.team_scratch(scratch_level));
      ViewTypeBScratch B_scr_real(team.team_scratch(scratch_level));
      ViewTypeBScratch B_scr_imag(team.team_scratch(scratch_level));
      ViewTypeCScratch C_scr(team.team_scratch(scratch_level));
      Kokkos::parallel_for(
         Kokkos::TeamThreadRange(team, blockA0), [&](const int i) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, blockB1),
               [&](const int j) { C_scr(i, j) = 0; });
         });
      team.team_barrier();

      // Move along the inner dimension in blocks
      const int length = TransposeA > 0 ? A.extent_int(0) : A.extent_int(1);
      for(int A_j = 0; A_j < length; A_j += blockA1)
      {
         // Load A block into scratch

         impl_deep_copy_matrix_block<
            typename Kokkos::TeamPolicy<ExecSpace>::member_type,
            ViewTypeAScratch, ViewTypeA,
            typename impl_gemm_choose_copy_layout<ExecSpace,
               typename ViewTypeA::array_layout,
               typename ViewTypeAScratch::array_layout>::type,
            blockA0, blockA1, TransposeA>::copy(team, A_scr, A, i_offset, A_j);

         // Load B block into scratch
         impl_deep_copy_matrix_block<
            typename Kokkos::TeamPolicy<ExecSpace>::member_type,
            ViewTypeBScratch, ViewTypeB,
            typename impl_gemm_choose_copy_layout<ExecSpace,
               typename ViewTypeB::array_layout,
               typename ViewTypeBScratch::array_layout>::type,
            blockA1, blockB1, TransposeB>::copy(team, B_scr_real, B_scr_imag, B, A_j, j_offset);

         // Wait for A and B block to be in scratch memory
         team.team_barrier();

         // Add contribution from multiplying the A and B block to the C block
         impl_team_gemm_block(team, C_scr, A_scr, B_scr_real);
         team.team_barrier();

         // Add contribution from multiplying the A and B block to the C block
         impl_team_gemm_block(team, C_scr, A_scr, B_scr_real);

         // Wait for subblock computation to be done before loading the next A
         // and B block
         team.team_barrier();
      }
      // Write back the C block from scratch to main memory
      impl_update_matrix_block<
         typename Kokkos::TeamPolicy<ExecSpace>::member_type, ViewTypeC,
         ViewTypeCScratch, typename ViewTypeC::array_layout, blockA0,
         blockB1>::update(team, beta, C, alpha, C_scr, i_offset, j_offset);
   }
};

template <typename ExecSpace> inline int get_max_vector_size() {
   return Kokkos::TeamPolicy<ExecSpace>::vector_length_max();
}

template <class AV, class BV, class CV>
struct GEMM {
   static void gemm(const typename CV::execution_space &space, const AV &A,
      const BV &B, const CV &C) {
      static_assert(
         Kokkos::is_view<AV>::value, "AV must be a Kokkos::View.");
      static_assert(
         Kokkos::is_view<BV>::value, "BV must be a Kokkos::View.");
      static_assert(
         Kokkos::is_view<CV>::value, "CV must be a Kokkos::View.");
      static_assert(
         static_cast<int>(AV::rank) == 2, "AV must have rank 2.");
      static_assert(
         static_cast<int>(BV::rank) == 2, "BV must have rank 2.");
      static_assert(
         static_cast<int>(CV::rank) == 2, "CV must have rank 2.");

      // Figure out Scalar Types
      typedef typename AV::non_const_value_type ScalarA;
      typedef typename BV::non_const_value_type ScalarB;
      typedef typename CV::non_const_value_type ScalarC;
      typedef typename CV::execution_space ExecSpace;

      // Define Blocking sizes (this will be used for scratch spaces)
      static constexpr int blockA0 = 24;
      static constexpr int blockB1 = 64;
      static constexpr int blockA1 =
         (sizeof(ScalarA) * blockA0 * 16 + sizeof(ScalarB) * 16 * blockB1 +
               sizeof(ScalarC) * blockA0 * blockB1 <
            24000)
            ? 16
         : (sizeof(ScalarA) * blockA0 * 8 + sizeof(ScalarB) * 8 * blockB1 +
                 sizeof(ScalarC) * blockA0 * blockB1 <
              24000)
            ? 8
         : (sizeof(ScalarA) * blockA0 * 4 + sizeof(ScalarB) * 4 * blockB1 +
                 sizeof(ScalarC) * blockA0 * blockB1 <
              24000)
            ? 4
            : 16;
      int vector_length = blockB1 / 4;
      int max_vector_length =
         get_max_vector_size<typename CV::execution_space>();
      if(vector_length > max_vector_length)
         vector_length = max_vector_length;

      // Compute scratch space size
      typedef GEMMImpl<typename CV::execution_space, AV, BV, CV, blockA0,
         blockA1, blockB1, 0, 0>
         gemm_dummy_type;
      const int scratch_memory_size =
         gemm_dummy_type::ViewTypeAScratch::required_allocation_size() +
         gemm_dummy_type::ViewTypeBScratch::required_allocation_size() +
         gemm_dummy_type::ViewTypeCScratch::required_allocation_size();
      const int scratch_level = scratch_memory_size < 24000 ? 0 : 1;

      // Figure out Team Sizes
      int team_size = 1;
#if defined(KOKKOS_ENABLE_CUDA)
      if(std::is_same<typename CV::execution_space, Kokkos::Cuda>::value)
         team_size = blockA0;
#endif
#if defined(KOKKOS_ENABLE_HIP)
      if(std::is_same<typename CV::execution_space,
            Kokkos::Experimental::HIP>::value)
         team_size = blockA0;
#endif
#if defined(KOKKOS_ENABLE_ROCM)
      if(std::is_same<typename CV::execution_space, Kokkos::ROCm>::value)
         team_size = blockA0;
#endif
#if defined(KOKKOS_ENABLE_SYCL)
      if(std::is_same<typename CV::execution_space,
            Kokkos::Experimental::SYCL>::value)
         team_size = blockA0;
#endif

      GEMMImpl<typename CV::execution_space, AV, BV, CV, blockA0, blockA1,
         blockB1, 0, 0>
         gemm(A, B, C);
      gemm.run(space, team_size, vector_length, scratch_level);
   }
};

template <class AViewType, class BViewType, class CViewType>
void gemm(const typename CViewType::execution_space& space, const char transA[],
          const char transB[], typename AViewType::const_value_type& alpha,
          const AViewType& A, const BViewType& B,
          typename CViewType::const_value_type& beta, const CViewType& C) {
  // Return if C matrix is degenerated
  if ((C.extent(0) == 0) || (C.extent(1) == 0)) {
    return;
  }

  // Minimize the number of Impl::GEMM instantiations, by
  // standardizing on particular View specializations for its template
  // parameters.
  typedef Kokkos::View<
      typename AViewType::const_value_type**, typename AViewType::array_layout,
      typename AViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      AVT;
  typedef Kokkos::View<
      typename BViewType::const_value_type**, typename BViewType::array_layout,
      typename BViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      BVT;
  typedef Kokkos::View<typename CViewType::non_const_value_type**,
                       typename CViewType::array_layout,
                       typename CViewType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      CVT;
  typedef GEMM<AVT, BVT, CVT> impl_type;
  impl_type::gemm(space, A, B, C);
}

template <class AViewType, class BViewType, class CViewType>
void gemm(const char transA[], const char transB[],
          typename AViewType::const_value_type& alpha, const AViewType& A,
          const BViewType& B, typename CViewType::const_value_type& beta,
          const CViewType& C) {
  const typename CViewType::execution_space space =
      typename CViewType::execution_space();
  gemm(space, transA, transB, alpha, A, B, beta, C);
}

namespace KokkosUnitOperator{

   //Operator working on horizontal rout layout
   template<typename R, typename T, typename V>
   struct ApplyOperator {

      ApplyOperator(R rv, R iv, T vo, V s, int tl, int ir, int cs)
          : rOutView(rv), inView(iv), Ops(vo), scan(s), total(tl), inRows(ir),
            col_size(cs) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int thid) const {
         auto local_row = thid % total;
         auto index = binary_search_range(scan, local_row);

         auto row = local_row - scan(index);
         auto col = thid / total + index * col_size;

         for(int j = 0; j < inRows; j++)
         { rOutView(row, col) += Ops(local_row, j) * inView(j, col); }
      }


      R rOutView;
      R inView;
      T Ops;
      V scan;
      int total;
      int inRows;
      int col_size;
   };


   //This works well rout and in vertical left but mops vertical right.
   template <typename R, typename L, typename T, typename V>
   struct ApplyUnitOperatorLeft {

      ApplyUnitOperatorLeft(R rv, L iv, T vo, V s, int ir)
          : rOutView(rv), inView(iv), Ops(vo), scan(s), inRows(ir) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int row, const int col) const {
         auto index = binary_search_range(scan, row);

         for(int j = 0; j < inRows; j++)
         {
            rOutView(row, col) += Ops(row, j) * inView(j + index * inRows, col);
         }
      }

      R rOutView;
      L inView;
      T Ops;
      V scan;
      int inRows;
   };


   template <typename R, typename L, typename T, typename V>
   struct ApplyUnitOperatorVerticalIN {

      ApplyUnitOperatorVerticalIN(R rv, L iv, T vo, V s, int ir)
          : rOutView(rv), inView(iv), Ops(vo), scan(s), inRows(ir) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int row, const int col) const {
         auto index = binary_search_range(scan, row);

         for(int j = 0; j < inRows; j++)
         {
            rOutView(row, col) += Ops(row, j) * inView(j + index * inRows, col);
         }
      }

      R rOutView;
      L inView;
      T Ops;
      V scan;
      int inRows;
   };


   template <typename R, typename L, typename T> struct ApplyUnitOperatorGemm {

      ApplyUnitOperatorGemm(R rv, L iv, T vo, int ir)
          : rOutView(rv), inView(iv), Ops(vo), inRows(ir) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int row, const int col) const {

         for(int j = 0; j < inRows; j++)
         { rOutView(row, col) += Ops(row, j) * inView(j, col); }
      }

      R rOutView;
      L inView;
      T Ops;
      int inRows;
   };


   //Requires rOut to be a vertical matrix instead of horizontal.
   //Better coalescing achieved with very good occupancy.
   //However requires specialized copy of rout into its original format
   template <typename R, typename L, typename T, typename V>
   struct ApplyUnitOperatorLoopExchange{

      ApplyUnitOperatorLoopExchange(R rv, L iv, T vo, V s, int cs)
          : rOutView(rv), inView(iv), Ops(vo), scan(s), col_size(cs) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int row, const int col) const {
         auto index = binary_search_range(scan, row);

         for(int j = 0; j < col_size; j++)
         {
            auto colin = j + index * col_size;
            rOutView(row, j) += Ops(row, col) * inView(col, colin);
         }
      }

      R rOutView;
      L inView;
      T Ops;
      V scan;
      int col_size;
   };

   template <typename R, typename L, typename T, typename V>
   struct ApplyUnitOperatorLoopExchangeVerticalIN{

      ApplyUnitOperatorLoopExchangeVerticalIN(R rv, L iv, T vo, V s, int cs, int inr)
          : rOutView(rv), inView(iv), Ops(vo), scan(s), col_size(cs), inRows(inr) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int row, const int col) const {
         auto index = binary_search_range(scan, row);

         for(int j = 0; j < col_size; j++)
         {
            auto colin = j + index * col_size;
            rOutView(row, j) += Ops(row, col) * inView(col + index * inRows, j);
         }
      }

      R rOutView;
      L inView;
      T Ops;
      V scan;
      int col_size;
      int inRows;
   };

   //Requires rOut to be a vertical matrix instead of horizontal.
   //Better coalescing achieved with very good occupancy.
   //However requires specialized copy of rout into its original format
   template <typename R, typename L, typename T, typename V>
   struct ApplyUnitOperatorLoopExchange1D{

      ApplyUnitOperatorLoopExchange1D(R rv, L iv, T vo, V s, int ir, int cs)
          : rOutView(rv), inView(iv), Ops(vo), scan(s), inRows(ir),
            col_size(cs) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int row) const {
         auto index = binary_search_range(scan, row);

         for(int k = 0; k < col_size; k++)
         {
            for(int j = 0; j < inRows; j++)
            {
               auto colin = j + index * col_size;
               rOutView(row, j) += Ops(row, k) * inView(k, colin);
            }
         }
      }

      R rOutView;
      L inView;
      T Ops;
      V scan;
      int inRows;
      int col_size;
   };


   template <typename T = double, typename CC = Kokkos::complex<double>>
   void test_gemm() {
      int M = 5000;
      int N = 4000;

      ViewMatrixType<T> A("A:", M, N);
      ViewMatrixTypeLeft<CC> B("B:", N, M);
      ViewMatrixType<T> C("C:", M, M);

      Kokkos::deep_copy(A, 1.0);
      Kokkos::deep_copy(B, 2.0);

      const T alpha = T(1.0);
      const T beta = T(0.0);

      KokkosBlas::gemm("N", "N", alpha, A, B, beta, C);

      gemm("N", "N", alpha, A, B, beta, C);

      Kokkos::parallel_for("Apply Kokkos Operator",
         Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
            {0, 0}, {C.extent(0), C.extent(1)}),
         ApplyUnitOperatorGemm(C, B, A, B.extent(0)));
   }


   template <typename OpTypes>
   void applyUnitOperatorTeam(const typename OpTypes::OpMatrixLZ &rOutView,
      const typename OpTypes::OpMatrixLZ &inView,
      const typename OpTypes::OpMatrixL Ops,
      const typename OpTypes::OpVectorI &scan, const int total,
      const int inRows, const int col_size) {

      using DataType = typename OpTypes::DataType;

      using member_type = Kokkos::TeamPolicy<>::member_type;
      using team_policy = Kokkos::TeamPolicy<>;

      //Same as the best performance without team policy
      //but needs layout right.
      //TODO: Check what happens with horizontal matrix for OPS.
      Kokkos::parallel_for(
         "gmm", team_policy(total, col_size),
         KOKKOS_LAMBDA(const member_type &teamMember) {
            const int local_col = teamMember.team_rank();
            const int local_row = teamMember.league_rank();
            DataType temp2 = 0;

            auto index = binary_search_range(scan, local_row);

            auto row = local_row - scan(index);
            auto col = local_col + index * col_size;

            for(int j = 0; j < inRows; j++)
            { rOutView(row, col) += Ops(local_row, j) * inView(j, col); }
         });

      //Reduction version which is not as performant as without.
      //Best performance so far is one thread per row.
      Kokkos::parallel_for(
         "gmm", team_policy(total * col_size, Kokkos::AUTO),
         "gmm", team_policy(total * col_size, 1),
         KOKKOS_LAMBDA(const member_type &teamMember) {
            const int thid = teamMember.team_rank();
            const int team = teamMember.league_rank();

            auto local_row = team % total;
            auto index = binary_search_range(scan, local_row);

            auto row = local_row - scan(index);
            auto col = team / total + index * col_size;

            for(int j = 0; j < inRows; j++)
            { rOutView(row, col) += Ops(local_row, j) * inView(j, col); }

            DataType temp2 = 0;
            Kokkos::parallel_reduce(
               Kokkos::TeamThreadRange(teamMember, inRows),
               [&](const int j, DataType &innerUpdate) {
                  innerUpdate += Ops(local_row, j) * inView(j, col);
               },
               temp2);

            if(thid == 0)
               rOutView(row, col) = temp2;
         });
   }

   template <typename OpMatrixLZ, typename OpMatrixLZL, typename OpMatrixL,
      typename OpVectorI>
   void testUnitOperator(const OpMatrixLZ& rOutView, const OpMatrixLZL& inView,
      const OpMatrixL& vmOps, const OpVectorI& scan, const int total,
      const int inRows, const int col_size)
   {
      int row_size = rOutView.extent(0);

      Kokkos::parallel_for("Apply1 Kokkos Operator",
         Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {row_size, inRows}),
         ApplyUnitOperatorLoopExchangeVerticalIN(rOutView, inView, vmOps, scan,
            col_size, inRows));

      Kokkos::parallel_for("Apply Kokkos Operator",
         Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0},
            {rOutView.extent(0), rOutView.extent(1)}),
         ApplyUnitOperatorVerticalIN(rOutView, inView, vmOps, scan, inRows));

      Kokkos::parallel_for("Apply Kokkos Operator",
         Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0},
            {rOutView.extent(0), rOutView.extent(1)}),
         ApplyUnitOperator(rOutView, inView, vmOps, scan, inRows, col_size));

      test_gemm();
   }
}
#endif
} // namespace Experimental
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_PIOPERATORTYPES_HPP
