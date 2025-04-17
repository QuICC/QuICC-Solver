/**
 * @file SphRadLapl.hpp
 * @brief Implementation of the Chebyshev based radial part of spherical
 * Laplacian 1/Y^2 D Y^2 D projector, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_SPHRADLAPL_BASE_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_SPHRADLAPL_BASE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/ILinearMapProjector.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Tags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

template <typename Impl> class SphRadLapl;

/**
 * @brief Implementation of the Chebyshev based radial part of spherical
 * Laplacian 1/Y^2 D Y^2 D projector, with linear map y = ax + b
 */
template <> class SphRadLapl<base_t> : public ILinearMapProjector
{
public:
   /**
    * @brief Constructor
    */
   SphRadLapl() = default;

   /**
    * @brief Destructor
    */
   virtual ~SphRadLapl() = default;

protected:
   /**
    * @brief Initialize storage
    */
   virtual void initBackend() const override;

private:
   /**
    * @brief Initialize solver operators
    */
   virtual void initOperator() const override;

   /**
    * @brief Apply pre FFT operator
    *
    * @param rOut Output values
    * @param in   Input values
    */
   virtual void applyPreOperator(Matrix& rOut, const Matrix& in) const override;

   /**
    * @brief Apply post FFT operator
    *
    * @param rOut Output values
    */
   virtual void applyPostOperator(Matrix& rOut) const override;

   /**
    * @brief Apply pre FFT operator for component wise openerations
    *
    * @param in   Input values
    * @param useReal Real vs Imag flag
    */
   virtual void applyPreOperator(Matrix& tmp, const MatrixZ& in,
      const bool useReal) const override;

   /**
    * @brief Apply post FFT operator for component wise operations
    *
    * @param rOut Output values
    * @param useReal Real vs Imag flag
    */
   virtual void applyPostOperator(MatrixZ& rOut, const Matrix& tmp,
      const bool useReal) const override;
};

} // namespace Projector
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_SPHRADLAPL_BASE_HPP
