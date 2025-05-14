/**
 * @file SphereTorPolEnstrophyLSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics
 * enstrophy L spectrum calculation for toroidal/poloidal field in a sphere
 */

// System includes
//

// Project includes
//
#include "QuICC/Io/Variable/SphereTorPolEnstrophyLSpectrumWriter.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

SphereTorPolEnstrophyLSpectrumWriter::SphereTorPolEnstrophyLSpectrumWriter(
   const std::string& prefix, const std::string& type) :
    ISphericalTorPolEnstrophyLSpectrumWriter(prefix, type)
{}

void SphereTorPolEnstrophyLSpectrumWriter::init()
{
   // Sphere volume: 4/3*pi*r_o^3
   this->mVolume = (4.0 / 3.0) * Math::PI;

   this->mHasMOrdering = this->res().sim().ss().has(
      SpatialScheme::Feature::TransformSpectralOrdering123);

   ISphericalTorPolEnstrophyLSpectrumWriter::init();
}

} // namespace Variable
} // namespace Io
} // namespace QuICC
