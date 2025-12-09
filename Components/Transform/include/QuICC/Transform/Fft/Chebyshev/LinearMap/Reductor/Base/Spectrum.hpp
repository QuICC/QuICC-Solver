/**
 * @file Spectrum.hpp
 * @brief Implementation of the Chebyshev based 1/Y^1 radial power reductor,
 * with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_SPECTRUM_BASE_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_SPECTRUM_BASE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/ILinearMapSpectrum.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Tags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

template <typename Impl> class Spectrum;

/**
 * @brief Implementation of the Chebyshev based 1/Y^1 radial power reductor,
 * with linear map y = ax + b
 */
template <> class Spectrum<base_t> : public ILinearMapSpectrum
{
public:
   /**
    * @brief Constructor
    */
   Spectrum() = default;

   /**
    * @brief Destructor
    */
   ~Spectrum() = default;

protected:
   /**
    * @brief Initialise operator
    */
   void initOperator() const final;

private:
   /**
    * @brief Apply pre FFT operator
    *
    * @param in   Input values
    */
   void applyPreOperator(Matrix& tmp, const Matrix& in) const final;

   /**
    * @brief Apply post FFT operator
    *
    * @param rOut Output values
    */
   void applyPostOperator(Matrix& rOut, const Matrix& tmp) const final;

   /**
    * @brief Apply pre FFT operator for component wise openerations
    *
    * @param in   Input values
    * @param useReal Real vs Imag flag
    */
   void applyPreOperator(Matrix& tmp, const MatrixZ& in,
      const bool useReal) const final;
};

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_SPECTRUM_BASE_HPP
