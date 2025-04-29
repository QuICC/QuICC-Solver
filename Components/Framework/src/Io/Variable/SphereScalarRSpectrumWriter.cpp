/**
 * @file SphereScalarRSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics radial
 * power spectrum calculation for scalar field in a sphere
 */

// System includes
//
#include <iomanip>

// Project includes
//
#include "QuICC/Io/Variable/SphereScalarRSpectrumWriter.hpp"
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/SphereScalarRSpectrumWriter.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

SphereScalarRSpectrumWriter::SphereScalarRSpectrumWriter(
   const std::string& prefix, const std::string& type) :
    ISphericalScalarRSpectrumWriter(prefix, type)
{}

void SphereScalarRSpectrumWriter::init()
{
   // Normalize by sphere volume: 4/3*pi*r_o^3
   this->mVolume = (4.0 / 3.0) * Math::PI;

   this->mHasMOrdering = this->res().sim().ss().has(
      SpatialScheme::Feature::TransformSpectralOrdering123);

   ISphericalScalarRSpectrumWriter::init();
}

} // namespace Variable
} // namespace Io
} // namespace QuICC
