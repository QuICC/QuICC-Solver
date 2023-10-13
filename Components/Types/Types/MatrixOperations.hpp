/**
 * @file MatrixOperations.hpp
 * @brief Useful methods for the DecoupledComplex type
 */

#ifndef QUICC_TYPES_MATRIXOPERATIONS_HPP
#define QUICC_TYPES_MATRIXOPERATIONS_HPP

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Datatypes {

namespace details {
//
// Matrix products
//

template <typename T1, typename T2, typename T3>
void setMatrixProduct(T1& rField, const int start, const T2& mat,
   const T3& rhs);

template <typename T1, typename T2, typename T3, typename T4>
void setMatrixProduct(T1& rField, const int start, const T2& mat,
   const T3& rhsRe, const T4& rhsIm);

template <typename T1, typename T2, typename T3>
void addMatrixProduct(T1& rField, const int start, const T2& mat,
   const T3& rhs);

template <typename T1, typename T2, typename T3, typename T4>
void addMatrixProduct(T1& rField, const int start, const T2& mat,
   const T3& rhsRe, const T4& rhsIm);

//
// 1D Top block operations
//

template <typename T1, typename T2>
void setTopBlock(T1& rField, const int start, const int rows, const T2& rhs);

//
// 2D Top block operations
//

template <typename T1, typename T2>
void setTopBlock(T1& rField, const int start, const int rows,
   const int fastSize, const int fastShift, const T2& rhs);

/*
 * Inline definitions below
 */

template <typename T1, typename T2, typename T3>
inline void setMatrixProduct(T1& rField, const int start, const T2& mat,
   const T3& rhs)
{
   if constexpr (std::is_same<T1, Matrix>::value ||
                 std::is_same<T1, MatrixZ>::value)
   {
      rField.block(start, 0, mat.rows(), rField.cols()) = mat * rhs;
   }
   else if constexpr (std::is_same<T1, DecoupledZMatrix>::value)
   {
      int rows = mat.rows();
      int cols = rField.real().cols();

      if constexpr (std::is_same<T2, SparseMatrix>::value)
      {
         rField.real().block(start, 0, rows, cols) = mat * rhs.real();
         rField.imag().block(start, 0, rows, cols) = mat * rhs.imag();
      }
      else if (std::is_same<T2, SparseMatrixZ>::value)
      {
         assert(rField.real().rows() == rField.imag().rows());
         assert(rField.real().cols() == rField.imag().cols());

         rField.real().block(start, 0, rows, cols) =
            mat.real() * rhs.real() - mat.imag() * rhs.imag();
         rField.imag().block(start, 0, rows, cols) =
            mat.real() * rhs.imag() + mat.imag() * rhs.real();
      }
      else
      {
         static_assert(true,
            "Tried to use invalid combination of types in setMatrixProduct");
      }
   }
   else
   {
      static_assert(true,
         "Tried to use invalid combination of types in setMatrixProduct");
   }
}

template <typename T1, typename T2, typename T3, typename T4>
inline void setMatrixProduct(T1& rField, const int start, const T2& mat,
   const T3& rhsRe, const T4& rhsIm)
{
   if constexpr (std::is_same<T1, DecoupledZMatrix>::value)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      int rows = mat.rows();
      int cols = rField.real().cols();
      if constexpr (std::is_same<T2, SparseMatrix>::value)
      {
         rField.real().block(start, 0, rows, cols) = mat * rhsRe;
         rField.imag().block(start, 0, rows, cols) = mat * rhsIm;
      }
      else if constexpr (std::is_same<T2, SparseMatrixZ>::value)
      {
         rField.real().block(start, 0, rows, cols) =
            mat.real() * rhsRe - mat.imag() * rhsIm;
         rField.imag().block(start, 0, rows, cols) =
            mat.real() * rhsIm + mat.imag() * rhsRe;
      }
      else
      {
         static_assert(true, "Tried to use invalid combination of types in "
                             "setMatrixProduct for split RSH");
      }
   }
   else
   {
      static_assert(true, "Tried to use invalid combination of types in "
                          "setMatrixProduct for split RSH");
   }
}

template <typename T1, typename T2, typename T3>
inline void addMatrixProduct(T1& rField, const int start, const T2& mat,
   const T3& rhs)
{
   if constexpr (std::is_same<T1, Matrix>::value ||
                 std::is_same<T1, MatrixZ>::value)
   {
      rField.block(start, 0, mat.rows(), rField.cols()) += mat * rhs;
   }
   else if constexpr (std::is_same<T1, DecoupledZMatrix>::value)
   {
      int rows = mat.rows();
      int cols = rField.real().cols();

      if constexpr (std::is_same<T2, SparseMatrix>::value)
      {
         rField.real().block(start, 0, rows, cols) += mat * rhs.real();
         rField.imag().block(start, 0, rows, cols) += mat * rhs.imag();
      }
      else if (std::is_same<T2, SparseMatrixZ>::value)
      {
         assert(rField.real().rows() == rField.imag().rows());
         assert(rField.real().cols() == rField.imag().cols());

         rField.real().block(start, 0, rows, cols) +=
            mat.real() * rhs.real() - mat.imag() * rhs.imag();
         rField.imag().block(start, 0, rows, cols) +=
            mat.real() * rhs.imag() + mat.imag() * rhs.real();
      }
      else
      {
         static_assert(true,
            "Tried to use invalid combination of types in addMatrixProduct");
      }
   }
   else
   {
      static_assert(true,
         "Tried to use invalid combination of types in addMatrixProduct");
   }
}

template <typename T1, typename T2, typename T3, typename T4>
inline void addMatrixProduct(T1& rField, const int start, const T2& mat,
   const T3& rhsRe, const T4& rhsIm)
{
   if constexpr (std::is_same<T1, DecoupledZMatrix>::value)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      int rows = mat.rows();
      int cols = rField.real().cols();
      if constexpr (std::is_same<T2, SparseMatrix>::value)
      {
         rField.real().block(start, 0, rows, cols) += mat * rhsRe;
         rField.imag().block(start, 0, rows, cols) += mat * rhsIm;
      }
      else if (std::is_same<T2, SparseMatrixZ>::value)
      {
         rField.real().block(start, 0, rows, cols) +=
            mat.real() * rhsRe - mat.imag() * rhsIm;
         rField.imag().block(start, 0, rows, cols) +=
            mat.real() * rhsIm + mat.imag() * rhsRe;
      }
      else
      {
         static_assert(true, "Tried to use invalid combination of types in "
                             "addMatrixProduct for split RHS");
      }
   }
   else
   {
      static_assert(true, "Tried to use invalid combination of types in "
                          "addMatrixProduct for split RHS");
   }
}

template <typename T1, typename T2>
inline void setTopBlock(T1& rField, const int start, const int rows,
   const T2& rhs)
{
   if constexpr (std::is_same<T1, Matrix>::value ||
                 std::is_same<T1, MatrixZ>::value)
   {
      if constexpr (std::is_same<typename T1::Scalar,
                       typename T2::Scalar>::value)
      {
         int cols = rField.cols();
         rField.block(start, 0, rows, cols) = rhs.topRows(rows);
      }
      else
      {
         static_assert(true,
            "Tried to use invalid combination of types in setTopBlock");
      }
   }
   else if constexpr (std::is_same<T1, DecoupledZMatrix>::value &&
                      std::is_same<T1, T2>::value)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      int cols = rField.real().cols();
      rField.real().block(start, 0, rows, cols) = rhs.real().topRows(rows);
      rField.imag().block(start, 0, rows, cols) = rhs.imag().topRows(rows);
   }
   else
   {
      static_assert(true,
         "Tried to use invalid combination of types in setTopBlock");
   }
}

template <typename T1, typename T2>
inline void setTopBlock(T1& rField, const int start, const int rows,
   const int fastSize, const int fastShift, const T2& rhs)
{
   if constexpr (std::is_same<T1, Matrix>::value ||
                 std::is_same<T1, MatrixZ>::value)
   {
      if constexpr (std::is_same<typename T1::Scalar,
                       typename T2::Scalar>::value)
      {
         int cols = rField.cols();
         int galBlock = fastSize - fastShift;
         int nJ = rows / galBlock;
         assert(rows - nJ * galBlock == 0);
         assert(rhs.rows() >= nJ * fastSize);
         for (int j = 0; j < nJ; j++)
         {
            rField.block(j * galBlock + start, 0, galBlock, cols) =
               rhs.block(j * fastSize, 0, galBlock, cols);
         }
      }
      else
      {
         static_assert(true,
            "Tried to use invalid combination of types in setTopBlock");
      }
   }
   else if constexpr (std::is_same<T1, DecoupledZMatrix>::value &&
                      std::is_same<T1, T2>::value)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      int cols = rField.real().cols();
      int galBlock = fastSize - fastShift;
      int nJ = rows / galBlock;
      assert(rows - nJ * galBlock == 0);
      assert(rhs.real().rows() >= nJ * fastSize);
      for (int j = 0; j < nJ; j++)
      {
         rField.real().block(j * galBlock + start, 0, galBlock, cols) =
            rhs.real().block(j * fastSize, 0, galBlock, cols);
         rField.imag().block(j * galBlock + start, 0, galBlock, cols) =
            rhs.imag().block(j * fastSize, 0, galBlock, cols);
      }
   }
   else
   {
      static_assert(true,
         "Tried to use invalid combination of types in setTopBlock");
   }
}
} // namespace details
} // namespace Datatypes
} // namespace QuICC

#endif // QUICC_TYPES_MATRIXOPERATIONS_HPP
