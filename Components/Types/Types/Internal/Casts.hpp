/**
 * @file Casts.hpp
 * @brief Casts for internal structures to base types
 * computations
 */

#ifndef QUICC_TYPES_INTERNAL_CASTS_HPP
#define QUICC_TYPES_INTERNAL_CASTS_HPP

// System includes
//
#include <type_traits>

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {
namespace Internal {

/// @brief Cast the internal value to an external one
/// @param val Internal::MHDFloat
/// @return QuICC::MHDFloat
inline QuICC::MHDFloat cast(Internal::MHDFloat val)
{
#ifdef QUICC_MULTPRECISION
   return val.convert_to<QuICC::MHDFloat>();
#else
   // flush to zero subnormal values
   if (std::abs(val) < std::numeric_limits<QuICC::MHDFloat>::min())
   {
      val = 0.0;
   }
   return val;
#endif // QUICC_MULTPRECISION
}

namespace details {

/// @brief generic Eigen Dense structure cast from internal Scalar to base
/// Scalar
template <class T1, class T2,
   std::enable_if_t<std::is_same_v<typename T1::Scalar, QuICC::MHDFloat> &&
                       std::is_same_v<typename T2::Scalar, Internal::MHDFloat>,
      bool> = true>
inline void cast(Eigen::DenseBase<T1>& OMat, const Eigen::DenseBase<T2>& IMat)
{
   // Loop over whole matrix
   for (int j = 0; j < IMat.cols(); ++j)
   {
      for (int i = 0; i < IMat.rows(); ++i)
      {
         OMat(i, j) = Internal::cast(IMat(i, j));
      }
   }
}
} // namespace details

/// @brief Cast the internal smart Array to an external one
/// @param spIArr Internal smart Array to cast
/// @return QuICC::SharedArray
inline QuICC::SharedArray cast(Internal::SharedArray spIArr)
{
   auto spArr = std::make_shared<QuICC::Array>(spIArr->size());

   details::cast(*spArr, *spIArr);

   return spArr;
}

/// @brief Cast the internal smart Matrix to an external one
/// @param spIMat Internal smart Matrix to cast
/// @return QuICC::SharedMatrix
inline QuICC::SharedMatrix cast(Internal::SharedMatrix spIMat)
{
   auto spMat = std::make_shared<QuICC::Matrix>(spIMat->rows(), spIMat->cols());

   details::cast(*spMat, *spIMat);

   return spMat;
}

/// @brief Cast the internal Eigen::Ref Array to an external one
/// @param rIArr Internal Eigen::Ref Array to cast
/// @return QuICC::Array
inline QuICC::Array cast(const Internal::Array& rIArr)
{
   QuICC::Array arr(rIArr.size());

   details::cast(arr, rIArr);

   return arr;
}

/// @brief Cast the internal Eigen::Ref Matrix to an external one
/// @param rIMat Internal Eigen::Ref Matrix to cast
/// @return QuICC::Matrix
inline QuICC::Matrix cast(const Internal::Matrix& rIMat)
{
   QuICC::Matrix mat(rIMat.rows(), rIMat.cols());

   details::cast(mat, rIMat);

   return mat;
}

/// @brief Cast the internal Eigen::Ref Matrix to an external one
/// @param rIMat Internal Eigen::Ref Matrix to cast
/// @return QuICC::Matrix
inline QuICC::Matrix cast(const Eigen::Ref<const Internal::Matrix>& rIMat)
{
   QuICC::Matrix mat(rIMat.rows(), rIMat.cols());

   details::cast(mat, rIMat);

   return mat;
}

} // namespace Internal
} // namespace QuICC

#endif // QUICC_TYPES_INTERNAL_CASTS_HPP
