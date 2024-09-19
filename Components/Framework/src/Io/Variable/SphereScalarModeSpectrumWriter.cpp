/**
 * @file SphereScalarModeSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics mode
 * energy spectrum calculation for scalar field in a sphere
 */

// System includes
//
#include <iomanip>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/SphereScalarModeSpectrumWriter.hpp"
#include "Environment/QuICCEnv.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

SphereScalarModeSpectrumWriter::SphereScalarModeSpectrumWriter(
   const std::string& prefix, const std::string& type) :
    ISphericalScalarModeSpectrumWriter(prefix, type)
{}

void SphereScalarModeSpectrumWriter::init()
{
   // Sphere volume (r_o = 1): 4/3*pi*r_o^3
   this->mVolume = (4.0 / 3.0) * Math::PI;

   this->mHasMOrdering = this->res().sim().ss().has(
      SpatialScheme::Feature::TransformSpectralOrdering123);

   ISphericalScalarModeSpectrumWriter::init();
}

} // namespace Variable
} // namespace Io
} // namespace QuICC
