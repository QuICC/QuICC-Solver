/**
 * @file Tags.hpp
 * @brief Fourier operators backends
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_TAGS_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_TAGS_HPP

// System includes
//

// Project includes
//

namespace QuICC {
namespace Transform {
namespace Fft {
namespace Fourier {

    /// old implementaion tag
    struct base_t {};

    /// view cpu implementation tag
    struct viewCpu_t {};

    /// view gpu implementation tag
    struct viewGpu_t {};

} // namespace Fourier
} // namespace Fft
} // namespace Transform
} // namespace QuICC


#endif // QUICC_TRANSFORM_FFT_FOURIER_TAGS_HPP
