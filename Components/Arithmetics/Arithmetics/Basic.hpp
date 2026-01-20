/**
 * @file Basic.hpp
 * @brief Useful generic methods for basis arithmetic operations on various types
 */

#ifndef QUICC_ARITHMETICS_BASIC_HPP
#define QUICC_ARITHMETICS_BASIC_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Arithmetics/Utility.hpp"
#include "View/View.hpp"

namespace QuICC {

namespace Arithmetics {

/**
 * @brief Get stored value
 *
 * @tparam TData
 * @param mat
 * @param k
 */
template <typename TData> auto getScalar(const TData& mat, const int k);

/**
 * @brief Get stored value in matrix
 *
 * @tparam TData
 * @param mat
 * @param i
 * @param j
 */
template <typename TData>
auto getScalar(const TData& mat, const int i, const int j);

/**
 * @brief Assign value
 *
 * @tparam TData
 * @param mat
 * @param k
 * @param val
 */
template <Operation EQUALOP, typename T1, typename T2>
void assignScalar(T1& mat, const int k, const T2& val);

/**
 * @brief Assign value
 *
 * @tparam TData
 * @param mat
 * @param i
 * @param j
 * @param val
 */
template <Operation EQUALOP, typename T1, typename T2>
void assignScalar(T1& mat, const int i, const int j, const T2& val);

/**
 * @brief Zero value
 *
 * @tparam TData
 * @param mat
 * @param k
 * @param val
 */
template <typename T1>
void setZero(T1& mat, const int k);

/**
 * @brief Zero value
 *
 * @tparam TData
 * @param mat
 * @param i
 * @param j
 * @param val
 */
template <typename T1>
void setZero(T1& mat, const int i, const int j);

namespace details {
   template <typename T>
   inline T& setS(Eigen::SparseMatrix<T>& mat, const int k)
   {
      return mat.coeffRef(k, 0);
   }

   template <typename T>
   inline T& setS(Eigen::SparseMatrix<T>& mat, const int i, const int j)
   {
      return mat.coeffRef(i, j);
   }

   template <typename T>
   inline T& setS(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const int k)
   {
      return mat(k);
   }

   template <typename T>
   inline T& setS(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const int i, const int j)
   {
      return mat(i, j);
   }

   template <typename T1, typename T2>
   inline T1& setS(View::View<T1,T2>& mat, const int k)
   {
      return mat[k];
   }

   template <typename T1, typename T2>
   inline T1& setS(View::View<T1,T2>& mat, const int i, const int j)
   {
      return mat(i, j);
   }

   template <typename T>
   inline T getS(const MHDVariant val)
   {
      return std::get<T>(val);
   }

   template <typename T>
   inline T getS(const T val)
   {
      return val;
   }
}

template <typename TData> inline auto getScalar(const TData& mat, const int k)
{
   if constexpr (std::is_same<TData, DecoupledZMatrix>::value)
   {
      return MHDComplex(mat.real()(k), mat.imag()(k));
   }
   else if constexpr(is_view<TData>::value)
   {
      return mat[k];
   }
   else
   {
      return mat(k);
   }
}

template <typename TData>
inline auto getScalar(const TData& mat, const int i, const int j)
{
   if constexpr (std::is_same<TData, DecoupledZMatrix>::value)
   {
      return MHDComplex(mat.real()(i, j), mat.imag()(i, j));
   }
   else
   {
      return mat(i, j);
   }
}

template <Operation EQUALOP, typename T1, typename T2>
inline void assignScalar(T1& mat, const int k, const T2& val)
{
   if constexpr (std::is_same<T1, DecoupledZMatrix>::value)
   {
      if constexpr(EQUALOP == Operation::Set)
      {
         mat.real()(k) = details::getS<MHDComplex>(val).real();
         mat.imag()(k) = details::getS<MHDComplex>(val).imag();
      }
      else if constexpr(EQUALOP == Operation::Plus)
      {
         mat.real()(k) += details::getS<MHDComplex>(val).real();
         mat.imag()(k) += details::getS<MHDComplex>(val).imag();
      }
      else if constexpr(EQUALOP == Operation::Minus)
      {
         mat.real()(k) -= details::getS<MHDComplex>(val).real();
         mat.imag()(k) -= details::getS<MHDComplex>(val).imag();
      }
   }
   else
   {
      if constexpr(EQUALOP == Operation::Set)
      {
         details::setS(mat, k) = details::getS<typename GetScalarType<T1>::ScalarType>(val);
      }
      else if constexpr(EQUALOP == Operation::Plus)
      {
         details::setS(mat, k) += details::getS<typename GetScalarType<T1>::ScalarType>(val);
      }
      else if constexpr(EQUALOP == Operation::Minus)
      {
         details::setS(mat, k) -= details::getS<typename GetScalarType<T1>::ScalarType>(val);
      }
   }
}

template <Operation EQUALOP, typename T1, typename T2>
inline void assignScalar(T1& mat, const int i, const int j, const T2& val)
{
   if constexpr (std::is_same<T1, DecoupledZMatrix>::value)
   {
      if constexpr(EQUALOP == Operation::Set)
      {
         mat.real()(i, j) = details::getS<MHDComplex>(val).real();
         mat.imag()(i, j) = details::getS<MHDComplex>(val).imag();
      }
      else if constexpr(EQUALOP == Operation::Plus)
      {
         mat.real()(i, j) += details::getS<MHDComplex>(val).real();
         mat.imag()(i, j) += details::getS<MHDComplex>(val).imag();
      }
      else if constexpr(EQUALOP == Operation::Minus)
      {
         mat.real()(i, j) -= details::getS<MHDComplex>(val).real();
         mat.imag()(i, j) -= details::getS<MHDComplex>(val).imag();
      }
   }
   else
   {
      if constexpr(EQUALOP == Operation::Set)
      {
         details::setS(mat, i, j) = details::getS<typename GetScalarType<T1>::ScalarType>(val);
      }
      else if constexpr(EQUALOP == Operation::Plus)
      {
         details::setS(mat, i, j) += details::getS<typename GetScalarType<T1>::ScalarType>(val);
      }
      else if constexpr(EQUALOP == Operation::Minus)
      {
         details::setS(mat, i, j) -= details::getS<typename GetScalarType<T1>::ScalarType>(val);
      }
   }
}

template <typename T1>
inline void setZero(T1& mat, const int k)
{
   assignScalar<Operation::Set>(mat, k, typename GetScalarType<T1>::ScalarType(0.0));
}

template <typename T1>
inline void setZero(T1& mat, const int i, const int j)
{
   assignScalar<Operation::Set>(mat, i, j, typename GetScalarType<T1>::ScalarType(0.0));
}

} // namespace Arithmetics
} // namespace QuICC

#endif // QUICC_ARITHMETICS_BASIC_HPP
