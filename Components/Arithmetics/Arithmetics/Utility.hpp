/**
 * @file Utility.hpp
 * @brief Useful generic utility methods for operations on various types
 */

#ifndef QUICC_ARITHMETICS_UTILITY_HPP
#define QUICC_ARITHMETICS_UTILITY_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "View/View.hpp"

namespace QuICC {

namespace Arithmetics {

enum class Operation { Minus = -1, Set = 0, Plus = 1};

/**
 * @brief Identify type of scalar
 */
template <typename TData> struct GetScalarType
{
   typedef typename TData::Scalar ScalarType;
};

/**
 * @brief Identify type of scalar of DecoupledComplex
 */
template <> struct GetScalarType<DecoupledZMatrix>
{
   typedef MHDComplex ScalarType;
};

/**
 * @brief Identify type of scalar of View
 */
template <typename T1, typename T2> struct GetScalarType<View::View<T1, T2>>
{
   typedef T1 ScalarType;
};

/**
 * @brief Is view?
 */
template <typename TData> struct is_view: std::false_type {};

/**
 * @brief Is view?
 */
template <typename T1, typename T2> struct is_view<View::View<T1, T2>>: std::true_type {};

/**
 * @brief Get number of columns
 */
template <typename TData> auto getCols(const TData& mat);

template <typename T1, typename T2> auto getCols(const View::View<T1,T2>& mat);

/**
 * @brief Class for holding a temporary
 */
template <typename TData>
struct Temporary
{
   const TData* ptr;
   TData  data;
};

template <typename TData> auto getCols(const TData& mat)
{
   if constexpr(std::is_same_v<TData, DecoupledZMatrix>)
   {
      return mat.real().cols();
   }
   else
   {
      return mat.cols();
   }
}

template <typename T1, typename T2> auto getCols(const View::View<T1,T2>& mat)
{
   return mat.dims()[1];
}

/**
 * @brief Class for holding temporary in View format
 */
template <typename T1, typename T2>
struct Temporary<View::View<T1,T2>>
{
   using TData = View::View<T1,T2>;
   const TData* ptr;
   TData  data;
   std::vector<T1> storage;
};

} // namespace Arithmetics
} // namespace QuICC

#endif // QUICC_ARITHMETICS_UTILITY_HPP
