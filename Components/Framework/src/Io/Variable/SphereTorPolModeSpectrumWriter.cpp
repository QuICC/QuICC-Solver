/**
 * @file SphereTorPolModeSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy
 * calculation for toroidal/poloidal field in a sphere
 */

// System includes
//
#include <iomanip>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/SphereTorPolModeSpectrumWriter.hpp"
#include "Environment/QuICCEnv.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

SphereTorPolModeSpectrumWriter::SphereTorPolModeSpectrumWriter(
   const std::string& prefix, const std::string& type) :
    ISphericalTorPolModeSpectrumWriter(prefix, type)
{}

void SphereTorPolModeSpectrumWriter::init()
{
   // Sphere volume: 4/3*pi*r_o^3
   this->mVolume = (4.0 / 3.0) * Math::PI;

   this->mHasMOrdering = this->res().sim().ss().has(
      SpatialScheme::Feature::TransformSpectralOrdering123);

   ISphericalTorPolModeSpectrumWriter::init();
}

} // namespace Variable
} // namespace Io
} // namespace QuICC
