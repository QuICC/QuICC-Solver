/**
 * @file MatrixOperations.hpp
 * @brief Useful methods for the DecoupledComplex type
 */

#ifndef QUICC_ARITHMETICS_LINEARALGEBRA_HPP
#define QUICC_ARITHMETICS_LINEARALGEBRA_HPP

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Arithmetics/Utility.hpp"

namespace QuICC {

namespace Arithmetics {

namespace details {

   template <Operation EQUALOP, typename TMat>
   void computeAx(MatrixZ& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const Matrix& rhsRe, const Matrix& rhsIm);

   template <Operation EQUALOP, Operation SIGN>
   void computeAxBy(Matrix& rField, const std::tuple<int,int,int,int>& outBlk, const SparseMatrix& mat, const Matrix& x, const SparseMatrix& matB, const Matrix& y);

   template <Operation EQUALOP, Operation SIGN>
   void computeAxBy(Matrix& rField, const std::tuple<int,int,int,int>& outBlk, const SparseMatrixZ& mat, const Matrix& x, const Matrix& y);

   template <Operation EQUALOP, Operation SIGN>
   void computeAxBy(Matrix& rField, const std::tuple<int,int,int,int>& outBlk, const SparseMatrixZ& mat, const Matrix& x, const Matrix& y, const std::tuple<int,int,int,int>& inBlk);
}

//
// Matrix products
//
template <Operation EQUALOP, typename TOut, typename TMat, typename TIn>
void computeAx(Eigen::Matrix<TOut, Eigen::Dynamic, Eigen::Dynamic>& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const Eigen::Matrix<TIn, Eigen::Dynamic, Eigen::Dynamic>& rhs);

template <Operation EQUALOP, typename TOut, typename TMat, typename TIn>
void computeAx(Eigen::Matrix<TOut, Eigen::Dynamic, Eigen::Dynamic>& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const Eigen::Matrix<TIn, Eigen::Dynamic, Eigen::Dynamic>& rhs, const std::tuple<int,int,int,int>& inBlk);

template <Operation EQUALOP, typename TMat>
void computeAx(DecoupledZMatrix& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const DecoupledZMatrix& rhs);

template <Operation EQUALOP, typename TMat>
void computeAx(DecoupledZMatrix& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const DecoupledZMatrix& rhs, const std::tuple<int,int,int,int>& inBlk);

template <Operation EQUALOP, typename TMat>
void computeAx(DecoupledZMatrix& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const MatrixZ& rhs);

template <Operation EQUALOP, typename TMat>
void computeAx(MatrixZ& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const DecoupledZMatrix& rhs);

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

template <typename T1, typename T2>
void setTopBlock(T1& rField, const int start, const int rows,
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
         static_assert(false,
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
      static_assert(false,
         "Tried to use invalid combination of types in setTopBlock");
   }
}

template <typename T1, typename T2>
void setTopBlock(T1& rField, const int start, const int rows,
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
         static_assert(false,
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
      static_assert(false,
         "Tried to use invalid combination of types in setTopBlock");
   }
}

   template <Operation EQUALOP, typename TOut, typename TMat, typename TIn>
   void computeAx(Eigen::Matrix<TOut, Eigen::Dynamic, Eigen::Dynamic>& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const Eigen::Matrix<TIn, Eigen::Dynamic, Eigen::Dynamic>& rhs)
   {
      if constexpr(std::is_same_v<TIn, TOut>)
      {
         if constexpr(EQUALOP == Operation::Set)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) = mat * rhs;
         }
         else if constexpr(EQUALOP == Operation::Plus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) += mat * rhs;
         }
         else if constexpr(EQUALOP == Operation::Minus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) -= mat * rhs;
         }
      }
      else
      {
         throw std::logic_error("Incompatible types in computeAx call");
      }
   }

   template <Operation EQUALOP, typename TOut, typename TMat, typename TIn>
   void computeAx(Eigen::Matrix<TOut, Eigen::Dynamic, Eigen::Dynamic>& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const Eigen::Matrix<TIn, Eigen::Dynamic, Eigen::Dynamic>& rhs, const std::tuple<int,int,int,int>& inBlk)
   {
      if constexpr(std::is_same_v<TIn, TOut>)
      {
         if constexpr(EQUALOP == Operation::Set)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) = mat * rhs.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk));
         }
         else if constexpr(EQUALOP == Operation::Plus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) += mat * rhs.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk));
         }
         else if constexpr(EQUALOP == Operation::Minus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) -= mat * rhs.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk));
         }
      }
      else
      {
         throw std::logic_error("Incompatible types in computeAx call (2)");
      }
   }

   template <Operation EQUALOP, typename TMat>
   void computeAx(DecoupledZMatrix& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const DecoupledZMatrix& rhs)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      auto& rFieldRe = rField.real();
      auto& rFieldIm = rField.imag();
      const auto& rhsRe = rhs.real();
      const auto& rhsIm = rhs.imag();
      if constexpr (std::is_same<TMat, SparseMatrix>::value)
      {
         computeAx<EQUALOP>(rFieldRe, outBlk, mat, rhsRe);
         computeAx<EQUALOP>(rFieldIm, outBlk, mat, rhsIm);
      }
      else if constexpr (std::is_same<TMat, SparseMatrixZ>::value)
      {
         details::computeAxBy<EQUALOP, Operation::Minus>(rFieldRe, outBlk, mat, rhsRe, rhsIm);
         details::computeAxBy<EQUALOP, Operation::Plus>(rFieldIm, outBlk, mat, rhsIm, rhsRe);
      }
      else
      {
         static_assert(false, "Tried to use invalid combination of types in "
               "setMatrixProduct for split RHS");
      }
   }

   template <Operation EQUALOP, typename TMat>
   void computeAx(DecoupledZMatrix& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const DecoupledZMatrix& rhs, const std::tuple<int,int,int,int>& inBlk)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      const auto& rhsRe = rhs.real();
      const auto& rhsIm = rhs.imag();
      if constexpr (std::is_same<TMat, SparseMatrix>::value)
      {
         computeAx<EQUALOP>(rField.real(), outBlk, mat, rhsRe, inBlk);
         computeAx<EQUALOP>(rField.imag(), outBlk, mat, rhsIm, inBlk);
      }
      else if constexpr (std::is_same<TMat, SparseMatrixZ>::value)
      {
         details::computeAxBy<EQUALOP, Operation::Minus>(rField.real(), outBlk, mat, rhsRe, rhsIm, inBlk);
         details::computeAxBy<EQUALOP, Operation::Plus>(rField.imag(), outBlk, mat, rhsIm, rhsRe, inBlk);
      }
      else
      {
         static_assert(false, "Tried to use invalid combination of types in "
               "setMatrixProduct for split RHS");
      }
   }

   template <Operation EQUALOP, typename TMat>
   void computeAx(DecoupledZMatrix& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const MatrixZ& rhs)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      auto& rFieldRe = rField.real();
      auto& rFieldIm = rField.imag();
      if constexpr (std::is_same<TMat, SparseMatrix>::value)
      {
         if constexpr(EQUALOP == Operation::Set)
         {
            rFieldRe.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) = mat * rhs.real();
            rFieldIm.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) = mat * rhs.imag();
         }
         else if constexpr(EQUALOP == Operation::Plus)
         {
            rFieldRe.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) += mat * rhs.real();
            rFieldIm.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) += mat * rhs.imag();
         }
         else if constexpr(EQUALOP == Operation::Minus)
         {
            rFieldRe.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) -= mat * rhs.real();
            rFieldIm.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) -= mat * rhs.imag();
         }
      }
      else if constexpr (std::is_same<TMat, SparseMatrixZ>::value)
      {
         if constexpr(EQUALOP == Operation::Set)
         {
            rFieldRe.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) = mat.real() * rhs.real() - mat.imag() * rhs.imag();
            rFieldIm.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) = mat.real() * rhs.imag() + mat.imag() * rhs.real();
         }
         else if constexpr(EQUALOP == Operation::Plus)
         {
            rFieldRe.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) += mat.real() * rhs.real() - mat.imag() * rhs.imag();
            rFieldIm.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) += mat.real() * rhs.imag() + mat.imag() * rhs.real();
         }
         else if constexpr(EQUALOP == Operation::Minus)
         {
            rFieldRe.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) -= mat.real() * rhs.real() - mat.imag() * rhs.imag();
            rFieldIm.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) -= mat.real() * rhs.imag() + mat.imag() * rhs.real();
         }
      }
      else
      {
         static_assert(false, "Tried to use invalid combination of types in "
               "setMatrixProduct for split RHS");
      }
   }

   template <Operation EQUALOP, typename TMat>
   void computeAx(MatrixZ& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const DecoupledZMatrix& rhs)
   {
      computeAx(rField, outBlk, mat, rhs.real(), rhs.imag());
   }

namespace details {
   template <Operation EQUALOP, typename TMat>
   void computeAx(MatrixZ& rField, const std::tuple<int,int,int,int>& outBlk, const TMat& mat, const Matrix& rhsRe, const Matrix& rhsIm)
   {
      if constexpr (std::is_same<TMat, SparseMatrix>::value)
      {
         if constexpr(EQUALOP == Operation::Set)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)).real() = mat * rhsRe;
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)).imag() = mat * rhsIm;
         }
         else if constexpr(EQUALOP == Operation::Plus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)).real() += mat * rhsRe;
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)).imag() += mat * rhsIm;
         }
         else if constexpr(EQUALOP == Operation::Minus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)).real() -= mat * rhsRe;
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)).imag() -= mat * rhsIm;
         }
      }
      else if constexpr (std::is_same<TMat, SparseMatrixZ>::value)
      {
         computeAxBy<EQUALOP, Operation::Minus>(rField.real(), outBlk, mat.real(), rhsRe, mat.imag(), rhsIm);
         computeAxBy<EQUALOP, Operation::Plus>(rField.imag(), outBlk, mat.real(), rhsIm, mat.imag(), rhsRe);
      }
      else
      {
         static_assert(false, "Tried to use invalid combination of types in "
               "setMatrixProduct for split RHS");
      }
   }

   template <Operation EQUALOP, Operation SIGN>
   void computeAxBy(Matrix& rField, const std::tuple<int,int,int,int>& outBlk, const SparseMatrix& matA, const Matrix& x, const SparseMatrix& matB, const Matrix& y)
   {
      assert(matA.rows() == matB.rows());
      assert(matA.cols() == matB.cols());

      if constexpr(SIGN == Operation::Plus)
      {
         if constexpr(EQUALOP == Operation::Set)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) = matA * x + matB * y;
         }
         else if constexpr(EQUALOP == Operation::Plus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) += matA * x + matB * y;
         }
         else if constexpr(EQUALOP == Operation::Minus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) -= matA * x + matB * y;
         }
      }
      else if constexpr(SIGN == Operation::Minus)
      {
         if constexpr(EQUALOP == Operation::Set)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) = matA * x - matB * y;
         }
         else if constexpr(EQUALOP == Operation::Plus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) += matA * x - matB * y;
         }
         else if constexpr(EQUALOP == Operation::Minus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) -= matA * x - matB * y;
         }
      }
      else
      {
         static_assert(false, "Undefined sign template parameter");
      }
   }

   template <Operation EQUALOP, Operation SIGN>
   void computeAxBy(Matrix& rField, const std::tuple<int,int,int,int>& outBlk, const SparseMatrixZ& mat, const Matrix& x, const Matrix& y)
   {
      if constexpr(SIGN == Operation::Plus)
      {
         if constexpr(EQUALOP == Operation::Set)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) = mat.real() * x + mat.imag() * y;
         }
         else if constexpr(EQUALOP == Operation::Plus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) += mat.real() * x + mat.imag() * y;
         }
         else if constexpr(EQUALOP == Operation::Minus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) -= mat.real() * x + mat.imag() * y;
         }
      }
      else if constexpr(SIGN == Operation::Minus)
      {
         if constexpr(EQUALOP == Operation::Set)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) = mat.real() * x - mat.imag() * y;
         }
         else if constexpr(EQUALOP == Operation::Plus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) += mat.real() * x - mat.imag() * y;
         }
         else if constexpr(EQUALOP == Operation::Minus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) -= mat.real() * x - mat.imag() * y;
         }
      }
      else
      {
         static_assert(false, "Undefined sign template parameter");
      }
   }

   template <Operation EQUALOP, Operation SIGN>
   void computeAxBy(Matrix& rField, const std::tuple<int,int,int,int>& outBlk, const SparseMatrixZ& mat, const Matrix& x, const Matrix& y, const std::tuple<int,int,int,int>& inBlk)
   {
      if constexpr(SIGN == Operation::Plus)
      {
         if constexpr(EQUALOP == Operation::Set)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) = mat.real() * x.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk)) + mat.imag() * y.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk));
         }
         else if constexpr(EQUALOP == Operation::Plus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) += mat.real() * x.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk)) + mat.imag() * y.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk));
         }
         else if constexpr(EQUALOP == Operation::Minus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) -= mat.real() * x.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk)) + mat.imag() * y.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk));
         }
      }
      else if constexpr(SIGN == Operation::Minus)
      {
         if constexpr(EQUALOP == Operation::Set)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) = mat.real() * x.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk)) - mat.imag() * y.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk));
         }
         else if constexpr(EQUALOP == Operation::Plus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) += mat.real() * x.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk)) - mat.imag() * y.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk));
         }
         else if constexpr(EQUALOP == Operation::Minus)
         {
            rField.block(std::get<0>(outBlk), std::get<1>(outBlk), std::get<2>(outBlk), std::get<3>(outBlk)) -= mat.real() * x.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk)) - mat.imag() * y.block(std::get<0>(inBlk), std::get<1>(inBlk), std::get<2>(inBlk), std::get<3>(inBlk));
         }
      }
      else
      {
         static_assert(false, "Undefined sign template parameter");
      }
   }
}

} // namespace Arithmetics
} // namespace QuICC

#endif // QUICC_ARITHMETICS_LINEARALGEBRA_HPP
