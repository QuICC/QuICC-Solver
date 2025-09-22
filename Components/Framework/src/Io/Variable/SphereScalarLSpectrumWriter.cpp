/**
 * @file SphereScalarLSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics L energy
 * spectrum calculation for scalar field in a sphere
 */

// System includes
//
#include <iomanip>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/SphereScalarLSpectrumWriter.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

SphereScalarLSpectrumWriter::SphereScalarLSpectrumWriter(
   const std::string& prefix, const std::string& type) :
    ISphericalScalarLSpectrumWriter(prefix, type)
{}

void SphereScalarLSpectrumWriter::init()
{
   // Normalize by sphere volume: 4/3*pi*r_o^3
   this->mVolume = (4.0 / 3.0) * Math::PI;

   this->mHasMOrdering = this->res().sim().ss().has(
      SpatialScheme::Feature::TransformSpectralOrdering123);

   ISphericalScalarLSpectrumWriter::init();
}

} // namespace Variable
} // namespace Io
} // namespace QuICC
