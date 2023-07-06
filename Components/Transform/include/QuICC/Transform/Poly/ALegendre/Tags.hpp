/**
 * @file Tags.hpp
 * @brief Associated Legendre operators backends
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_TAGS_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_TAGS_HPP

// System includes
//

// Project includes
//

namespace QuICC {
namespace Transform {
namespace Poly {
namespace ALegendre {

    /// old implementaion tag
    struct base_t {};

    /// kokkos API (kokkos and Cuda implementation) tag
    struct kokkos_t {};

    /// view cpu wrapper tag
    struct viewCpu_t {};

    /// view gpu wrapper tag
    struct viewGpu_t {};

} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC


#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_TAGS_HPP
