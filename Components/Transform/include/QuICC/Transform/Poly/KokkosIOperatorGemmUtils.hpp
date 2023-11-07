/**
 * @file KokkosIOperatorGemmUtils.hpp
 * @brief Associated Poly based operator kokkos irregular block gemm utils
 */

#ifndef QUICC_TRANSFORM_POLY_KOKKOSIOPERATORGEMMUTILS_HPP
#define QUICC_TRANSFORM_POLY_KOKKOSIOPERATORGEMMUTILS_HPP

// System includes
//

// Project includes
//

#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/KokkosCudaIOperatorGemmUtils.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

/* using DataType = cuDoubleComplex; */

// GEMM UTILS
//
// Thread block size
#define ROW_BLOCK_SIZE 16
#define COL_BLOCK_SIZE 16

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
   static void copy(const TeamHandle &team, const ViewTypeScratch &A_scr_real,
      const ViewTypeScratch &A_scr_imag, const ViewType &A, const int &offset_i,
      const int &offset_j, const int row_block_end, const int col_block_end) {
      if(offset_i + blockDim_i <= row_block_end &&
         offset_j + blockDim_j <= col_block_end)
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
                     auto complex_value = A(idx_i, idx_j);
                     A_scr_real(i, j) = complex_value.real();
                     A_scr_imag(i, j) = complex_value.imag();
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
                     auto complex_value =
                        idx_i < row_block_end && idx_j < col_block_end
                           ? A(idx_i, idx_j)
                           : ATV::zero();
                     A_scr_real(i, j) = complex_value.real();
                     A_scr_imag(i, j) = complex_value.imag();
                  });
            });
      }
   }

   KOKKOS_INLINE_FUNCTION
   static void copy(const TeamHandle &team, const ViewTypeScratch &A_scr,
      const ViewType &A, const int &offset_i, const int &offset_j,
      const int row_block_end, const int col_block_end) {
      if(offset_i + blockDim_i <= row_block_end &&
         offset_j + blockDim_j <= col_block_end)
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
                        idx_i < row_block_end && idx_j < col_block_end
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
   static void copy(const TeamHandle &team, const ViewTypeScratch &A_scr,
      const ViewType &A, const int &offset_i, const int &offset_j,
      const int row_block_end, const int col_block_end) {
      if(offset_i + blockDim_i <= row_block_end &&
         offset_j + blockDim_j <= col_block_end)
      {
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, blockDim_i), [&](const int i) {
               const int idx_i = offset_i + i;
               Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, blockDim_j),
                  [&](const int j) {
                     const int idx_j = offset_j + j;
                     A_scr(i, j) = A(idx_i, idx_j);
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
                     A_scr(i, j) =
                        idx_i < row_block_end && idx_j < col_block_end
                           ? A(idx_i, idx_j)
                           : ATV::zero();
                  });
            });
      }
   }

   KOKKOS_INLINE_FUNCTION
   static void copy(const TeamHandle &team, const ViewTypeScratch &A_scr_real,
      const ViewTypeScratch &A_scr_imag, const ViewType &A, const int &offset_i,
      const int &offset_j, const int row_block_end, const int col_block_end) {
      if(offset_i + blockDim_i <= row_block_end &&
         offset_j + blockDim_j <= col_block_end)
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
                        idx_i < row_block_end && idx_j < col_block_end
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
      const ViewType &A, const value_type &alpha,
      const ViewTypeScratch &A_scr_real, const ViewTypeScratch &A_scr_imag,
      const int &offset_i, const int &offset_j, const int row_block_end,
      const int col_block_end) {
      if(offset_i + blockDim_i <= row_block_end &&
         offset_j + blockDim_j <= col_block_end)
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
                        A(idx_i, idx_j) =
                           value_type(A_scr_real(i, j), A_scr_imag(i, j));
                     });
               } else
               {
                  Kokkos::parallel_for(
                     Kokkos::ThreadVectorRange(team, blockDim_i),
                     [&](const int i) {
                        const int idx_i = offset_i + i;
                        A(idx_i, idx_j) =
                           A(idx_i, idx_j) +
                           value_type(A_scr_real(i, j), A_scr_imag(i, j));
                     });
               }
            });
      } else
      {
         const int range_i = offset_i + blockDim_i <= row_block_end
                                ? blockDim_i
                                : row_block_end - offset_i;
         const int range_j = offset_j + blockDim_j <= col_block_end
                                ? blockDim_j
                                : col_block_end - offset_j;
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, range_j), [&](const int j) {
               const int idx_j = offset_j + j;
               if(beta == ATV::zero())
               {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, range_i),
                     [&](const int i) {
                        const int idx_i = offset_i + i;
                        A(idx_i, idx_j) =
                           value_type(A_scr_real(i, j), A_scr_imag(i, j));
                     });
               } else
               {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, range_i),
                     [&](const int i) {
                        const int idx_i = offset_i + i;
                        A(idx_i, idx_j) =
                           A(idx_i, idx_j) +
                           value_type(A_scr_real(i, j), A_scr_imag(i, j));
                     });
               }
            });
      }
   }
   KOKKOS_INLINE_FUNCTION
   static void update(const TeamHandle &team, const value_type &beta,
      const ViewType &A, const value_type &alpha, const ViewTypeScratch &A_scr,
      const int &offset_i, const int &offset_j, const int row_block_end,
      const int col_block_end) {
      if(offset_i + blockDim_i <= row_block_end &&
         offset_j + blockDim_j <= col_block_end)
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
                        A(idx_i, idx_j) = A_scr(i, j);
                     });
               } else
               {
                  Kokkos::parallel_for(
                     Kokkos::ThreadVectorRange(team, blockDim_i),
                     [&](const int i) {
                        const int idx_i = offset_i + i;
                        A(idx_i, idx_j) =
                           A(idx_i, idx_j) + A_scr(i, j);
                     });
               }
            });
      } else
      {
         const int range_i = offset_i + blockDim_i <= row_block_end
                                ? blockDim_i
                                : row_block_end - offset_i;
         const int range_j = offset_j + blockDim_j <= col_block_end
                                ? blockDim_j
                                : col_block_end - offset_j;
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, range_j), [&](const int j) {
               const int idx_j = offset_j + j;
               if(beta == ATV::zero())
               {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, range_i),
                     [&](const int i) {
                        const int idx_i = offset_i + i;
                        A(idx_i, idx_j) = A_scr(i, j);
                     });
               } else
               {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, range_i),
                     [&](const int i) {
                        const int idx_i = offset_i + i;
                        A(idx_i, idx_j) =
                           A(idx_i, idx_j) + A_scr(i, j);
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
      const ViewType &A, const value_type &alpha,
      const ViewTypeScratch &A_scr_real, const ViewTypeScratch &A_scr_imag,
      const int &offset_i, const int &offset_j, const int row_block_end,
      const int col_block_end) {
      if(offset_i + blockDim_i <= row_block_end &&
         offset_j + blockDim_j <= col_block_end)
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
                        A(idx_i, idx_j) =
                           value_type(A_scr_real(i, j), A_scr_imag(i, j));
                     });
               } else
               {
                  Kokkos::parallel_for(
                     Kokkos::ThreadVectorRange(team, blockDim_j),
                     [&](const int j) {
                        const int idx_j = offset_j + j;
                        A(idx_i, idx_j) =
                           A(idx_i, idx_j) + value_type(A_scr_real(i, j), A_scr_imag(i, j));
                     });
               }
            });
      } else
      {
         const int range_i = offset_i + blockDim_i <= row_block_end
                                ? blockDim_i
                                : row_block_end - offset_i;
         const int range_j = offset_j + blockDim_j <= col_block_end
                                ? blockDim_j
                                : col_block_end - offset_j;
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, range_i), [&](const int i) {
               const int idx_i = offset_i + i;
               if(beta == ATV::zero())
               {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, range_j),
                     [&](const int j) {
                        const int idx_j = offset_j + j;
                        A(idx_i, idx_j) = value_type(A_scr_real(i, j), A_scr_imag(i, j));
                     });
               } else
               {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, range_j),
                     [&](const int j) {
                        const int idx_j = offset_j + j;
                        A(idx_i, idx_j) =
                           A(idx_i, idx_j) +
                           value_type(A_scr_real(i, j), A_scr_imag(i, j));
                     });
               }
            });
      }
   }

   KOKKOS_INLINE_FUNCTION
   static void update(const TeamHandle &team, const value_type &beta,
      const ViewType &A, const value_type &alpha, const ViewTypeScratch &A_scr,
      const int &offset_i, const int &offset_j, const int row_block_end,
      const int col_block_end) {
      if(offset_i + blockDim_i <= row_block_end &&
         offset_j + blockDim_j <= col_block_end)
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
                        A(idx_i, idx_j) = A_scr(i, j);
                     });
               } else
               {
                  Kokkos::parallel_for(
                     Kokkos::ThreadVectorRange(team, blockDim_j),
                     [&](const int j) {
                        const int idx_j = offset_j + j;
                        A(idx_i, idx_j) =
                           A(idx_i, idx_j) + A_scr(i, j);
                     });
               }
            });
      } else
      {
         const int range_i = offset_i + blockDim_i <= row_block_end
                                ? blockDim_i
                                : row_block_end - offset_i;
         const int range_j = offset_j + blockDim_j <= col_block_end
                                ? blockDim_j
                                : col_block_end - offset_j;
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, range_i), [&](const int i) {
               const int idx_i = offset_i + i;
               if(beta == ATV::zero())
               {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, range_j),
                     [&](const int j) {
                        const int idx_j = offset_j + j;
                        A(idx_i, idx_j) = A_scr(i, j);
                     });
               } else
               {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, range_j),
                     [&](const int j) {
                        const int idx_j = offset_j + j;
                        A(idx_i, idx_j) = A(idx_i, idx_j) + A_scr(i, j);
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
   const ViewTypeC &C, const ViewTypeA &A, const ViewTypeB &BR,
   const ViewTypeB &BI) {
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
                  C_ij0 += A_ij * ScalarC(BR(j, B_j), BI(j, B_j));
                  C_ij1 += A_ij * ScalarC(BR(j, B_j + blockB1 / 4),
                                     BI(j, B_j + blockB1 / 4));
                  C_ij2 += A_ij * ScalarC(BR(j, B_j + 2 * blockB1 / 4),
                                     BI(j, B_j + 2 * blockB1 / 4));
                  C_ij3 += A_ij * ScalarC(BR(j, B_j + 3 * blockB1 / 4),
                                     BI(j, B_j + 3 * blockB1 / 4));
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
   class ViewS, int blockA0, int blockA1, int blockB1, int TransposeA,
   int TransposeB>
struct GEMMImpl {
   ViewTypeA A;
   ViewTypeB B;
   ViewTypeC C;
   ViewS allScan;
   ViewS scan;
   ViewS xGrid;
   ViewS yGrid;

   typedef typename ViewTypeA::non_const_value_type ScalarA;
   typedef typename ViewTypeB::non_const_value_type ScalarB;
   typedef typename ViewTypeC::non_const_value_type ScalarC;

   int num_blocks;
   /* const int num_blocks_1; */
   int scratch_level;

   ScalarC alpha, beta;
   typedef Kokkos::View<ScalarA[blockA0][blockA1], Kokkos::LayoutLeft,
      typename ExecSpace::scratch_memory_space>
      ViewTypeAScratch;

   typedef Kokkos::View<
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      ScalarA[blockA1][blockB1], Kokkos::LayoutRight,
#else
      ScalarA[blockA1][blockB1 + 1], Kokkos::LayoutRight,
#endif
      typename ExecSpace::scratch_memory_space>
      ViewTypeBScratch;

   typedef Kokkos::View<
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
       ScalarC[blockA0][blockB1], Kokkos::LayoutRight,
#else
       ScalarA[blockA0][blockB1 + 1], Kokkos::LayoutRight,
#endif
      typename ExecSpace::scratch_memory_space>
      ViewTypeCScratch;

   GEMMImpl(const ViewTypeA &A_, const ViewTypeB &B_, const ViewTypeC &C_,
      const ViewS &AS_, const ViewS &S_, const ViewS &XG_, const ViewS &YG_,
      const int nb_)
       : A(A_), B(B_), C(C_), allScan(AS_), scan(S_), xGrid(XG_), yGrid(YG_),
         /* num_blocks_0((C.extent_int(0) + blockA0 - 1) / blockA0), */
         num_blocks(nb_) {
      scratch_level = 0;
      alpha = 1;
      beta = 1;
   }

   GEMMImpl(const ScalarA &alpha_, const ViewTypeA &A_, const ViewTypeB &B_,
      const ScalarC &beta_, const ViewTypeC &C_, const int nb_)
       : A(A_), B(B_), C(C_),
         /* num_blocks_0((C.extent_int(0) + blockA0 - 1) / blockA0), */
         num_blocks(nb_) {
      scratch_level = 0;
      alpha = alpha_;
      beta = beta_;
   }

   void run(
      const ExecSpace &space, int team_size, int vector_length, int scr_level) {
      scratch_level = scr_level;
      int scratch_memory_size = ViewTypeAScratch::shmem_size() +
                                2 * ViewTypeBScratch::shmem_size() +
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
                                ViewTypeCScratch::shmem_size();
#else
                                2 * ViewTypeCScratch::shmem_size();
#endif


#if defined(KOKKOS_ENABLE_HIP)
      // Note lbv, 10/29/20: The LaunchBounds<384,2> leads
      // to an error with HIP as the heuristics on that platform
      // yield an optimal_num_blocks=0 which means no ressources
      // are allocated... Switching to LaunchBounds<384,2> fixes
      // that problem but I'm not sure if that it a good perf
      // parameter or why it is set to 2 for Cuda?
      Kokkos::TeamPolicy<ExecSpace, Kokkos::LaunchBounds<384, 0>> policy(
         space, num_blocks, team_size, vector_length);
#else
      Kokkos::TeamPolicy<ExecSpace, Kokkos::LaunchBounds<256, 2>> policy(
         space, num_blocks, team_size, vector_length);
#endif

      Kokkos::Timer timer_gen;
      Kokkos::parallel_for(impl_gemm_label<TransposeA, TransposeB>::label,
         policy.set_scratch_size(
            scratch_level, Kokkos::PerTeam(scratch_memory_size)),
         *this);
      double time_gen = timer_gen.seconds();
      std::cout << "KOKKOS operator teampolicy took: " << time_gen << " seconds."
                << std::endl;
   }

   KOKKOS_INLINE_FUNCTION
   void operator()(
      const typename Kokkos::TeamPolicy<ExecSpace>::member_type &team) const {
      // This team is responsible for computing a single block of C
      const int blockId = team.league_rank();

      // The id of the matrix block within the large matrix
      auto matrix_block_id = binary_search_range(allScan, blockId);
      // The start address of the A & B matrices
      // 2D block coordinates of each block id
      int blockRow = xGrid[blockId];
      int blockCol = yGrid[blockId];
      auto matrix_block_row_start = scan(matrix_block_id) + blockRow * blockA0;
      auto matrix_block_row_end = scan(matrix_block_id + 1);

      auto matrix_block_col_start = blockB1 * blockCol;
      auto matrix_blockB_start = matrix_block_id * C.extent(1);
      auto matrix_block_col_end = (matrix_block_id + 1) * C.extent(1);

      ViewTypeAScratch A_scr(team.team_scratch(scratch_level));
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      /* ViewTypeBScratch B_scr(team.team_scratch(scratch_level)); */
      ViewTypeBScratch B_scr_imag(team.team_scratch(scratch_level));
      ViewTypeBScratch B_scr_real(team.team_scratch(scratch_level));
      ViewTypeCScratch C_scr(team.team_scratch(scratch_level));
#else
      ViewTypeBScratch B_scr_imag(team.team_scratch(scratch_level));
      ViewTypeBScratch B_scr_real(team.team_scratch(scratch_level));
      ViewTypeCScratch C_scr_real(team.team_scratch(scratch_level));
      ViewTypeCScratch C_scr_imag(team.team_scratch(scratch_level));
#endif
      Kokkos::parallel_for(
         Kokkos::TeamThreadRange(team, blockA0), [&](const int i) {
            Kokkos::parallel_for(
               Kokkos::ThreadVectorRange(team, blockB1), [&](const int j) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
                  C_scr(i, j) = 0;
#else
                  C_scr_real(i, j) = 0;
                 C_scr_imag(i, j) = 0;
#endif
               });
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
            blockA0, blockA1, TransposeA>::copy(team, A_scr, A,
            matrix_block_row_start, A_j, matrix_block_row_end, A.extent(1));

         // Load B block into scratch
         impl_deep_copy_matrix_block<
            typename Kokkos::TeamPolicy<ExecSpace>::member_type,
            ViewTypeBScratch, ViewTypeB,
            typename impl_gemm_choose_copy_layout<ExecSpace,
               typename ViewTypeB::array_layout,
               typename ViewTypeBScratch::array_layout>::type,
            blockA1, blockB1, TransposeB>::copy(team, B_scr_real, B_scr_imag, B,
            A_j, matrix_block_col_start + matrix_blockB_start, B.extent(0),
            matrix_block_col_end);

         // Wait for A and B block to be in scratch memory
         team.team_barrier();

         // Add contribution from multiplying the A and B block to the C block
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
         /* impl_team_gemm_block(team, C_scr, A_scr, B_scr); */
         impl_team_gemm_block_complex(
            team, C_scr, A_scr, B_scr_real, B_scr_imag);
#else
         impl_team_gemm_block(team, C_scr_real, A_scr, B_scr_real);
         impl_team_gemm_block(team, C_scr_imag, A_scr, B_scr_imag);
#endif
         // Wait for subblock computation to be done before loading the next A
         // and B block
         team.team_barrier();
      }
      // Write back the C block from scratch to main memory
      impl_update_matrix_block<
         typename Kokkos::TeamPolicy<ExecSpace>::member_type, ViewTypeC,
         ViewTypeCScratch, typename ViewTypeC::array_layout, blockA0,
         blockB1>::update(team, beta, C, alpha,
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
         C_scr,
#else
         C_scr_real, C_scr_imag,
#endif
         matrix_block_row_start, matrix_block_col_start, matrix_block_row_end,
         C.extent(1));
   }
};

template <typename ExecSpace> inline int get_max_vector_size() {
   return Kokkos::TeamPolicy<ExecSpace>::vector_length_max();
}

template <class AV, class BV, class CV, class V>
struct GEMM {
   static void gemm(const typename CV::execution_space &space, const AV &A,
      const BV &B, const CV &C, const V &allScan, const V &scan, const V &xGrid,
      const V &yGrid, const int num_blocks) {
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
      static constexpr int blockA0 = ROW_BLOCK_SIZE;
      static constexpr int blockB1 = COL_BLOCK_SIZE;
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
      typedef GEMMImpl<typename CV::execution_space, AV, BV, CV, V, blockA0,
         blockA1, blockB1, 0, 0>
         gemm_dummy_type;
      const int scratch_memory_size =
         gemm_dummy_type::ViewTypeAScratch::required_allocation_size() +
         2 * gemm_dummy_type::ViewTypeBScratch::required_allocation_size() +
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
         gemm_dummy_type::ViewTypeCScratch::required_allocation_size();
#else
         2 * gemm_dummy_type::ViewTypeCScratch::required_allocation_size();
#endif
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

      GEMMImpl<typename CV::execution_space, AV, BV, CV, V, blockA0, blockA1,
         blockB1, 0, 0>
         gemm(A, B, C, allScan, scan, xGrid, yGrid, num_blocks);
      gemm.run(space, team_size, vector_length, scratch_level);
   }
};


template <class AViewType, class BViewType, class CViewType, class V>
void blockGemm(const AViewType &A, const BViewType &B, const CViewType &C,
   const V &allScan, const V &scan, const V &xGrid, const V &yGrid,
   const int allTotal) {
   const typename CViewType::execution_space space =
      typename CViewType::execution_space();
   // Return if C matrix is degenerated
   if((C.extent(0) == 0) || (C.extent(1) == 0))
   { return; }

   // Minimize the number of Impl::GEMM instantiations, by
   // standardizing on particular View specializations for its template
   // parameters.
   typedef Kokkos::View<typename AViewType::const_value_type **,
      typename AViewType::array_layout, typename AViewType::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      AVT;
   typedef Kokkos::View<typename BViewType::const_value_type **,
      typename BViewType::array_layout, typename BViewType::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      BVT;
   typedef Kokkos::View<typename CViewType::non_const_value_type **,
      typename CViewType::array_layout, typename CViewType::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      CVT;
   typedef GEMM<AVT, BVT, CVT, V> impl_type;
   impl_type::gemm(space, A, B, C, allScan, scan, xGrid, yGrid, allTotal);
}

template <Integer S = 0, typename T, typename MS, typename MZ, typename MZL,
   typename V>
void applyKokkosBlockOperator(const T mspSetup, const MS &vmOps, const MZ &rOutView,
   const MZL &inView, const V &scan, const int total) {

   auto slowSize = mspSetup->slowSize();
   auto outRows = mspSetup->fwdSize();

   V rowScan("outRows Scan", slowSize + 1);
   V colScan("outRows Scan", slowSize + 1);
   V allScan("outRows Scan", slowSize + 1);
   auto hostRowScan = Kokkos::create_mirror_view(rowScan);
   auto hostColScan = Kokkos::create_mirror_view(colScan);
   auto hostAllScan = Kokkos::create_mirror_view(allScan);

   static constexpr int blockA0 = ROW_BLOCK_SIZE;
   static constexpr int blockB1 = COL_BLOCK_SIZE;
   // build row and cols scan for each matrix using the block sizes.
   for(int i = 0; i < slowSize; i++)
   {
      if constexpr(S != 1)
      {
         outRows = mspSetup->fastSize(i);
      }
      auto col_size = mspSetup->mult(i);

      auto ro = outRows % blockA0;
      auto rowBlocks = outRows / blockA0;
      auto rc = col_size % blockB1;
      auto colBlocks = col_size / blockB1;

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

   // call the
   Kokkos::Timer timer_gen;
   blockGemm(vmOps, inView, rOutView, allScan, scan, xGrid, yGrid, allTotal);
   /* iKokkosBlockGemm<S>(vmOps, inView, rOutView, rowScan, colScan, scan, xGrid, yGrid,
      allScan, allTotal); */
   double time_gen = timer_gen.seconds();
   std::cout << "Kokkos ApplyBlockOp took: " << time_gen << " seconds."
             << std::endl;
}

#endif
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PIALEGENDREOPERATORTYPES_HPP
