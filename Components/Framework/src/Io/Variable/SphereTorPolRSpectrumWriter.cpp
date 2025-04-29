/**
 * @file SphereTorPolRSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics radial
 * power calculation for toroidal/poloidal field in a sphere
 */

// System includes
//
#include <iomanip>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/SphereTorPolRSpectrumWriter.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

SphereTorPolRSpectrumWriter::SphereTorPolRSpectrumWriter(
   const std::string& prefix, const std::string& type) :
    ISphericalTorPolRSpectrumWriter(prefix, type)
{}

void SphereTorPolRSpectrumWriter::init()
{
   // Sphere volume: 4/3*pi*r_o^3
   this->mVolume = (4.0 / 3.0) * Math::PI;

   this->mHasMOrdering = this->res().sim().ss().has(
      SpatialScheme::Feature::TransformSpectralOrdering123);

   ISphericalTorPolRSpectrumWriter::init();
}

} // namespace Variable
} // namespace Io
} // namespace QuICC
