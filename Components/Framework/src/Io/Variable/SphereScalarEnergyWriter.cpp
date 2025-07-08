/**
 * @file SphereScalarEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy
 * calculation for scalar field in a sphere
 */

// System includes
//
#include <iomanip>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/SphereScalarEnergyWriter.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

SphereScalarEnergyWriter::SphereScalarEnergyWriter(const std::string& prefix,
   const std::string& type) :
    ISphericalScalarEnergyWriter(prefix, type)
{}

void SphereScalarEnergyWriter::init()
{
   // Normalize by sphere volume: 4/3*pi*r_o^3
   this->mVolume = (4.0 / 3.0) * Math::PI;

   this->mHasMOrdering = this->res().sim().ss().has(
      SpatialScheme::Feature::TransformSpectralOrdering123);

   ISphericalScalarEnergyWriter::init();
}

} // namespace Variable
} // namespace Io
} // namespace QuICC
