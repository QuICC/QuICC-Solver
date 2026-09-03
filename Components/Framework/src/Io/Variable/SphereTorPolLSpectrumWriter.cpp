/**
 * @file SphereTorPolLSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy
 * calculation for toroidal/poloidal field in a sphere
 */

// System includes
//
#include <iomanip>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/SphereTorPolLSpectrumWriter.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

SphereTorPolLSpectrumWriter::SphereTorPolLSpectrumWriter(
   const std::string& prefix, const std::string& type, std::vector<std::shared_ptr<QuICC::DenseSM::Worland::RadialTorPolFunction>> pF) :
    ISphericalTorPolLSpectrumWriter(prefix, type, std::vector<std::shared_ptr<QuICC::DenseSM::IGenericProfile>>(pF.begin(), pF.end()))
{}

void SphereTorPolLSpectrumWriter::init()
{
   // Sphere volume: 4/3*pi*r_o^3
   this->mVolume = (4.0 / 3.0) * Math::PI;

   this->mHasMOrdering = this->res().sim().ss().has(
      SpatialScheme::Feature::TransformSpectralOrdering123);

   ISphericalTorPolLSpectrumWriter::init();
}

} // namespace Variable
} // namespace Io
} // namespace QuICC
