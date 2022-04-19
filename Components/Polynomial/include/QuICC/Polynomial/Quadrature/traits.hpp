/**
 * @file traits.hpp
 * @brief Implementation of the traits for quadratures
 */

#ifndef QUICC_POLYNOMIAL_QUADRATURE_TRAITS_HPP
#define QUICC_POLYNOMIAL_QUADRATURE_TRAITS_HPP

#include <type_traits>

namespace QuICC {
namespace Polynomial {
namespace Quadrature {

// normalization tags
static constexpr struct natural_t {} natural{};
static constexpr struct unity_t {} Unity{};

// Check normalization tags
template <class T>
struct isNormalizationTag : std::integral_constant
    <
        bool,
        std::is_same<natural_t, typename std::remove_cv<T>::type>::value ||
        std::is_same<unity_t, typename std::remove_cv<T>::type>::value
    > {};

// Generic template handles types that have no nested ::normalizationTag_t type member:
template <class, class = void>
struct hasNormalizationTag_t : std::false_type { };

// Specialization recognizes types that do have a nested ::normalizationTag_t type member:
template <class T>
struct hasNormalizationTag_t<T, std::void_t<typename T::normalizationTag_t>> : std::true_type {};



// quadrature type tags
static constexpr struct gauss_t {} gauss{};
static constexpr struct gaussLobatto_t {} gaussLobatto{};

// Check point type tags
template <class T>
struct isQuadTypeTag : std::integral_constant
    <
        bool,
        std::is_same<gauss_t, typename std::remove_cv<T>::type>::value ||
        std::is_same<gaussLobatto_t, typename std::remove_cv<T>::type>::value
    > {};

// Generic template handles types that have no nested ::quadTypeTag_t type member:
template <class, class = void>
struct hasQuadTypeTag_t : std::false_type { };

// Specialization recognizes types that do have a nested ::quadTypeTag_t type member:
template <class T>
struct hasQuadTypeTag_t<T, std::void_t<typename T::quadTypeTag_t>> : std::true_type {};


} // namespace Quadrature
} // namespace Polynomial
} // namespace QuICC


#endif // QUICC_POLYNOMIAL_QUADRATURE_TRAITS_HPP
