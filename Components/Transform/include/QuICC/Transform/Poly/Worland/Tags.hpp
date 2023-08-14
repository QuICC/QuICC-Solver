/**
 * @file Tags.hpp
 * @brief Worland operators backends
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_TAGS_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_TAGS_HPP

// System includes
//

// Project includes
//

namespace QuICC {
namespace Transform {
namespace Poly {
namespace Worland {

    /// old implementaion tag
    struct base_t {};

    /// kokkos API (kokkos and Cuda implementation) tag
    struct kokkos_t {};

    /// view cpu wrapper tag
    struct viewCpu_t {};

    /// view gpu wrapper tag
    struct viewGpu_t {};

} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC


#endif // QUICC_TRANSFORM_POLY_WORLAND_TAGS_HPP
