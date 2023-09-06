/**
 * @file KokkosUtils.hpp
 * @brief Some general Kokkos Utils
 */

#ifndef QUICC_KOKKOS_UTILS_HPP
#define QUICC_KOKKOS_UTILS_HPP

// Configuration includes
//

// System includes
//

#include "KokkosTypedefs.hpp"
#include<KokkosBlas3_gemm.hpp>
#include "QuICC/Transform/Poly/ALegendre/KokkosIALegendreOperatorTypes.hpp"

namespace QuICC {

/* special implementation of the binary search considering found
if an element is between current and next proc value. */
template <typename V, typename T>
KOKKOS_INLINE_FUNCTION Integer binary_search_range(
   const V view, const T value) {
   const int last_index = view.extent(0) - 1;
   int first_proc = 0;
   int last_proc = last_index;
   auto guess = last_index / 2;

   while(first_proc <= last_proc && first_proc != last_index)
   {
      T current = view(guess);
      T next = view(guess + 1);

      if(value >= current && value < next)
      {
         return guess;
      } else if(value < current)
      {
         last_proc = guess - 1;
         guess = (first_proc + last_proc + 1) / 2;
      } else if(value >= next)
      {
         first_proc = guess + 1;
         guess = (first_proc + last_proc) / 2;
      }
   }

   return -1;
}

// standard binary search
template <typename T>
KOKKOS_INLINE_FUNCTION Integer binary_search(
   const T *view, Integer f, Integer l, const T value) {
   while(f <= l)
   {
      Integer guess = (l + f) / 2;

      T current = *(view + guess);

      if(value == current)
      {
         return guess;
      } else if(value < current)
      {
         l = guess - 1;
      } else
      { f = guess + 1; }
   }

   return -1;
}

template <typename V, typename E>
void DeepCopyEigen(const V &view, const Eigen::DenseBase<E> &matrix) {
   auto rows = matrix.rows();
   auto cols = matrix.cols();

   using dataType = typename E::Scalar;
   const dataType *data = matrix.derived().data();

   ViewMatrixTypeLeftU<const dataType> matrix_view(data, rows, cols);

   auto hostView = Kokkos::create_mirror_view(view);

   Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>, KokkosHostExecSpace>(
         {0, 0}, {rows, cols}),
      [=](int i, int j) { hostView(i, j) = matrix_view(i, j); });

   Kokkos::deep_copy(view, hostView);
}

//assuming col size remains the same through each eigen matrix block.
template <typename V, typename E>
void DeepCopyEigen(
   const V &view, const Eigen::DenseBase<E> &matrix, const int col_size) {
   auto rows = matrix.rows();
   auto cols = matrix.cols();

   using dataType = typename E::Scalar;
   const dataType *data = matrix.derived().data();

   ViewMatrixTypeLeftU<const dataType> matrix_view(data, rows, cols);

   auto hostView = Kokkos::create_mirror_view(view);

   Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>, KokkosHostExecSpace>(
         {0, 0}, {rows, cols}),
      [=](int i, int j) {
         auto block = j / col_size;
         auto col = j % col_size;
         hostView(rows * block + i, col) = matrix_view(i, j);
      });

   Kokkos::deep_copy(view, hostView);
}

template <typename E, typename V>
void DeepCopyEigen(
   Eigen::DenseBase<E> &matrix, const V &view, const int col_size) {
   auto vrows = matrix.rows();
   auto vcols = matrix.cols();

   using dataType = typename E::Scalar;
   dataType *data = matrix.derived().data();

   ViewMatrixTypeLeftU<dataType> complex_copy(data, vrows, vcols);

   auto hostView = Kokkos::create_mirror_view(view);
   Kokkos::deep_copy(hostView, view);

   Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>, KokkosHostExecSpace>(
         {0, 0}, {vrows, vcols}),
      [=](int local_row, int local_col) {
         auto block = local_col / col_size;
         auto col = local_col % col_size;
         auto row = vrows * block + local_row;
         complex_copy(local_row, local_col).real(hostView(row, col).real());
         complex_copy(local_row, local_col).imag(hostView(row, col).imag());
      });
}

template <typename V, typename E>
void DeepCopyEigen(
   Eigen::DenseBase<E> &matrix, const V &view) {
   auto rows = matrix.rows();
   auto cols = matrix.cols();

   assert(rows == view.extent(0));
   assert(cols == view.extent(1));

   using dataType = typename E::Scalar;
   dataType *data = matrix.derived().data();

   ViewMatrixTypeLeftU<dataType> complex_copy(data, rows, cols);

   auto hostView = Kokkos::create_mirror_view(view);
   Kokkos::deep_copy(hostView, view);

   // Do a copy in parallel using OpenMP on the host views
   Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>, KokkosHostExecSpace>(
         {0, 0}, {rows, cols}),
      [=](int i, int j) { complex_copy(i, j) = hostView(i, j); });
}


template <typename VW, typename E, typename V, typename D>
void DeepCopyEigen(Eigen::DenseBase<E> &matrix,
   const VW &view, const V &scan, const D col_size) {
   using dataType = typename E::Scalar;
   dataType *data = matrix.derived().data();
   ViewMatrixTypeLeftU<dataType> complex_copy(data, matrix.rows(), matrix.cols());

   //This only works if the "view" has the same data type as eigen matrix (like dc).
   //meaning std::complex. If the view has kokkos::complex data it does not
   //Useful if the layout of the rout is switched to vertical. Then the rest is not needed.
   /* ViewMatrixTypeLeft<dataType> dc("DC", matrix.rows(), matrix.cols());
   Kokkos::deep_copy(complex_copy, view); */

   auto hostView = Kokkos::create_mirror_view(view);
   Kokkos::deep_copy(hostView, view);

   auto vrows = view.extent(0);
   auto vcols = view.extent(1);
   // Do a copy in parallel using OpenMP on the host views from a vertical
   // view to a horizontal (default Quicc layout) matrix.
   Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>, KokkosHostExecSpace>(
         {0, 0}, {vrows, vcols}),
      [=](int local_row, int local_col) {
         auto index = binary_search_range(scan, local_row);
         auto row = local_row - scan(index);
         auto col = index * col_size + local_col;
         complex_copy(row, col).real(hostView(local_row, local_col).real());
         complex_copy(row, col).imag(hostView(local_row, local_col).imag());
      });
}

template <typename T, typename E, typename V, typename D>
void DeepCopyEigen(T &view, const Eigen::DenseBase<E> &matrix, const V &scan,
   const D col_size) {
   using dataType = typename E::Scalar;
   const dataType *data = matrix.derived().data();
   ViewMatrixTypeLeftU<const dataType> complex_copy(data, matrix.rows(), matrix.cols());

   //This only works if the "view" has the same data type as eigen matrix (like dc).
   //meaning std::complex. If the view has kokkos::complex data it does not
   //Useful if the layout of the rout is switched to vertical. Then the rest is not needed.
   /* ViewMatrixTypeLeft<dataType> dc("DC", matrix.rows(), matrix.cols());
   Kokkos::deep_copy(complex_copy, view); */

   auto hostView = Kokkos::create_mirror_view(view);

   auto vrows = view.extent(0);
   auto vcols = view.extent(1);
   // Do a copy in parallel using OpenMP on the host views from a vertical
   // view to a horizontal (default Quicc layout) matrix.
   Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>, KokkosHostExecSpace>(
         {0, 0}, {vrows, vcols}),
      [=](int local_row, int local_col) {
         auto index = binary_search_range(scan, local_row);
         auto row = local_row - scan(index);
         auto col = index * col_size + local_col;
         hostView(local_row, local_col).real(complex_copy(row, col).real());
         hostView(local_row, local_col).imag(complex_copy(row, col).imag());
      });

   Kokkos::deep_copy(view, hostView);
}

template <typename S, typename T, typename E>
void denseMatrixMultiply(const S &SA, const S &SB,
   const ViewMatrixTypeLeft<T> &A, const ViewMatrixTypeLeft<E> &B,
   const ViewMatrixTypeLeft<E> &C) {

   const T alpha = T(1.0);
   const E beta = E(0.0);

   KokkosBlas::gemm(SA, SB, alpha, A, B, beta, C);
}


template <typename V, typename M>
void add_contribution_to_view_left(
   const V &hostView, const int index, const M &matrix) {

   Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>, KokkosHostExecSpace>(
         {0, 0}, {matrix.rows(), matrix.cols()}),
      [=](int i, int j) {
         auto l = i + index;
         hostView(l, j) = matrix(i, j);
      });
}

template <typename V, typename M>
void add_contribution_to_view_right(
   const V &hostView, const int index, const M &matrix) {

   Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>, KokkosHostExecSpace>(
         {0, 0}, {matrix.rows(), matrix.cols()}),
      [=](int i, int j) {
         auto l = j + index;
         hostView(i, l) = matrix(i, j);
      });
}


template <typename T, typename C>
bool complex_equal(T a, C b) {
   return (std::abs(a.real() - b.real()) < DBL_EPSILON &&
           std::abs(a.imag() - b.imag()) < DBL_EPSILON);
}

/* template <typename T>
bool complex_equal(T a, T b) {
   [>return std::abs(a - b) < DBL_EPSILON;<]
   return std::abs(a - b) < std::numeric_limits<double>::epsilon();
} */

/* bool complex_equal(double a, double b) {
   return (std::abs(a - b) <= DBL_EPSILON * std::max(std::abs(a), std::abs(b)));
} */

/* bool nearly_equal(double a, double b)
{
  return std::nextafter(a, std::numeric_limits<double>::lowest()) <= b
    && std::nextafter(a, std::numeric_limits<double>::max()) >= b;
} */

//check if matrix and view containt the same data. If not print all the differences.
template <typename M, typename V>
bool equal_data(const M &matrix, const V &view) {
   bool eq = true;
   for(int l = 0; l < view.extent(0); l++)
   {
      for(int k = 0; k < view.extent(1); k++)
      {
         if(!complex_equal(matrix(l, k), view(l, k)))
         /* if(!complex_equal(matrix(l, k).real(), view(l, k).real()) ||
            !complex_equal(matrix(l, k).imag(), view(l, k).imag())) */
         {
            eq = false;
            /* std::cout << "index: " << l << " " << k << " matrix: ("
                      << matrix(l, k).real() << " - " << matrix(l, k).imag()
                      << ") view (" << view(l, k).real() << " - "
                      << view(l, k).imag() << ")" << std::endl; */
         }
      }
   }
   return eq;
}


/*
bool equal_data(const ViewMatrixTypeLeft<MHDComplex> &matrix,
   const ViewMatrixTypeLeft<MHDComplex> &view) {

   for(int l = 0; l < matrix.extent(0); l++)
   {
      for(int k = 0; k < matrix.extent(1); k++)
      {
         if(!complex_equal(matrix(l, k), view(l, k)))
         {
             return false;
         }
      }
   }
   return true;
} */

template <typename T, typename V, typename S, typename D>
bool equal_data(
   const T &matrix, const V &view, const S &OutRowScan, D cols) {

   auto slowSize = OutRowScan.extent(0) - 1;

   bool eq = true;
   for(int l = 0; l < slowSize; l++)
   {
      auto outRows = OutRowScan(l + 1) - OutRowScan(l);
      for(int row = 0; row < outRows; row++)
      {
         for(int k = 0; k < cols; k++)
         {
            auto col = l * cols + k;

            if(!complex_equal(matrix(row, col), view(row, col)))
            {
               eq = false;
               /* std::cout << "index: " << row << " " << col << " matrix: ("
                         << matrix(row, col).real() << " - "
                         << matrix(row, col).imag() << ") view ("
                         << view(row, col).real() << " - "
                         << view(row, col).imag() << ")" << std::endl; */
            }
         }
      }
   }
   return eq;
}

} // namespace QuICC

#endif // QUICC_TYPEDEFS_HPP
