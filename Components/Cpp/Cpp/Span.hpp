/**
 * @file Span.hpp
 * @brief span wrapper for pre c++20
 */
#pragma once

// External includes
//
#ifdef __cpp_lib_span
#include <span>
namespace QuICC {
namespace Patch {
namespace std {

template <class T, ::std::size_t Extent = ::std::dynamic_extent>
    using span = ::std::span<T, Extent>;

} // namespace std
} // namespace Patch
} // namespace QuICC
#else
#include <boost/core/span.hpp>
namespace QuICC {
namespace Patch {
namespace std {

/// Typedef for boost implementation of span (until c++20)
template <class T, ::std::size_t E = boost::dynamic_extent>
    using span = boost::span<T, E>;

} // namespace std
} // namespace Patch
} // namespace QuICC
#endif
