/**
 * @file Tags.hpp
 * @brief Chebyshev LinearMap operators backends
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_TAGS_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_TAGS_HPP

// System includes
//

// Project includes
//

namespace QuICC {
namespace Transform {
namespace Fft {
namespace Chebyshev {
namespace LinearMap {

    /// old implementaion tag
    struct base_t {};

    /// view cpu wrapper tag
    struct viewCpu_t {};

    /// view gpu wrapper tag
    struct viewGpu_t {};

} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC


#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_TAGS_HPP
