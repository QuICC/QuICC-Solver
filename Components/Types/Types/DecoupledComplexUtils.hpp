/**
 * @file DecoupledComplexUtils.hpp
 * @brief Useful methods for the DecoupledComplex type
 */

#ifndef QUICC_TYPES_DECOUPLEDCOMPLEX_UTILS_HPP
#define QUICC_TYPES_DECOUPLEDCOMPLEX_UTILS_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Datatypes {

namespace details {
template <typename TData> class GetScalar
{
public:
   typedef typename TData::Scalar Scalar;
};

template <> class GetScalar<DecoupledZMatrix>
{
public:
   typedef MHDComplex Scalar;
};

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
 * @brief Set value
 *
 * @tparam TData
 * @param mat
 * @param k
 * @param val
 */
template <typename T1, typename T2>
void setScalar(T1& mat, const int k, const T2& val);

/**
 * @brief Set value
 *
 * @tparam TData
 * @param mat
 * @param i
 * @param j
 * @param val
 */
template <typename T1, typename T2>
void setScalar(T1& mat, const int i, const int j, const T2& val);

/**
 * @brief Add value
 *
 * @param mat
 * @param k
 * @param val
 */
template <typename T1, typename T2>
void addScalar(T1& mat, const int k, const T2& val);

/**
 * @brief Add value to DecoupledComplex storage
 *
 * @param mat
 * @param i
 * @param j
 * @param val
 */
template <typename T1, typename T2>
void addScalar(T1& mat, const int i, const int j, const T2& val);

template <typename TData> inline auto getScalar(const TData& mat, const int k)
{
   if constexpr (std::is_same<TData, DecoupledZMatrix>::value)
   {
      return MHDComplex(mat.real()(k), mat.imag()(k));
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

template <typename T1, typename T2>
inline void setScalar(T1& mat, const int k, const T2& val)
{
   if constexpr (std::is_same<T1, DecoupledZMatrix>::value &&
                 std::is_same<T2, MHDComplex>::value)
   {
      mat.real()(k) = val.real();
      mat.imag()(k) = val.imag();
   }
   else if constexpr (std::is_same<T1, DecoupledZMatrix>::value &&
                      std::is_same<T2, MHDVariant>::value)
   {
      mat.real()(k) = std::get<1>(val).real();
      mat.imag()(k) = std::get<1>(val).imag();
   }
   else if constexpr (std::is_same<T2, MHDVariant>::value &&
                      (std::is_same<T1,
                          Eigen::SparseMatrix<std::variant_alternative<0,
                             MHDVariant>::type>>::value ||
                         std::is_same<T1,
                            Eigen::SparseMatrix<std::variant_alternative<1,
                               MHDVariant>::type>>::value))
   {
      mat.coeffRef(k, 0) = std::get<typename T1::Scalar>(val);
   }
   else if constexpr (std::is_same<T1, Eigen::SparseMatrix<T2>>::value)
   {
      mat.coeffRef(k, 0) = val;
   }
   else if constexpr (std::is_same<T2, MHDVariant>::value)
   {
      mat(k) = std::get<typename T1::Scalar>(val);
   }
   else if constexpr (std::is_same<T2, typename T1::Scalar>::value)
   {
      mat(k) = val;
   }
   else
   {
      static_assert(true,
         "Tried to use invalid combination of types in setScalar(k)");
   }
}

template <typename T1, typename T2>
inline void setScalar(T1& mat, const int i, const int j, const T2& val)
{
   if constexpr (std::is_same<T1, DecoupledZMatrix>::value &&
                 std::is_same<T2, MHDComplex>::value)
   {
      mat.real()(i, j) = val.real();
      mat.imag()(i, j) = val.imag();
   }
   else if constexpr (std::is_same<T1, DecoupledZMatrix>::value &&
                      std::is_same<T2, MHDVariant>::value)
   {
      mat.real()(i, j) = std::get<1>(val).real();
      mat.imag()(i, j) = std::get<1>(val).imag();
   }
   else if constexpr (std::is_same<T2, MHDVariant>::value &&
                      (std::is_same<T1,
                          Eigen::SparseMatrix<std::variant_alternative<0,
                             MHDVariant>::type>>::value ||
                         std::is_same<T1,
                            Eigen::SparseMatrix<std::variant_alternative<1,
                               MHDVariant>::type>>::value))
   {
      mat.coeffRef(i, j) = std::get<typename T1::Scalar>(val);
   }
   else if constexpr (std::is_same<T1, Eigen::SparseMatrix<T2>>::value)
   {
      mat.coeffRef(i, j) = val;
   }
   else if constexpr (std::is_same<T2, MHDVariant>::value)
   {
      mat(i, j) = std::get<typename T1::Scalar>(val);
   }
   else if constexpr (std::is_same<T2, typename T1::Scalar>::value)
   {
      mat(i, j) = val;
   }
   else
   {
      static_assert(true,
         "Tried to use invalid combination of types in setScalar(i,j)");
   }
}

template <typename T1, typename T2>
inline void addScalar(T1& mat, const int k, const T2& val)
{
   if constexpr (std::is_same<T1, DecoupledZMatrix>::value &&
                 std::is_same<T2, MHDComplex>::value)
   {
      mat.real()(k) += val.real();
      mat.imag()(k) += val.imag();
   }
   else if constexpr (std::is_same<T1, DecoupledZMatrix>::value &&
                      std::is_same<T2, MHDVariant>::value)
   {
      mat.real()(k) += std::get<1>(val).real();
      mat.imag()(k) += std::get<1>(val).imag();
   }
   else if constexpr (std::is_same<T2, MHDVariant>::value)
   {
      mat(k) += std::get<typename T1::Scalar>(val);
   }
   else if constexpr (std::is_same<T2, typename T1::Scalar>::value)
   {
      mat(k) += val;
   }
   else
   {
      static_assert(true,
         "Tried to use invalid combination of types in addScalar(k)");
   }
}

template <typename T1, typename T2>
inline void addScalar(T1& mat, const int i, const int j, const T2& val)
{
   if constexpr (std::is_same<T1, DecoupledZMatrix>::value &&
                 std::is_same<T2, MHDComplex>::value)
   {
      mat.real()(i, j) += val.real();
      mat.imag()(i, j) += val.imag();
   }
   else if constexpr (std::is_same<T1, DecoupledZMatrix>::value &&
                      std::is_same<T2, MHDVariant>::value)
   {
      mat.real()(i, j) += std::get<1>(val).real();
      mat.imag()(i, j) += std::get<1>(val).imag();
   }
   else if constexpr (std::is_same<T2, MHDVariant>::value)
   {
      mat(i, j) += std::get<typename T1::Scalar>(val);
   }
   else if constexpr (std::is_same<T2, typename T1::Scalar>::value)
   {
      mat(i, j) += val;
   }
   else
   {
      static_assert(true,
         "Tried to use invalid combination of types in addScalar(i,j)");
   }
}
} // namespace details
} // namespace Datatypes
} // namespace QuICC

#endif // QUICC_TYPES_DECOUPLEDCOMPLEX_UTILS_HPP
