/**
 * @file SphereTorPolNCoeffsPowerWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics N spectra calculation for Tor and Pol scalar fields in a sphere
 */

// Configuration includes
//

// System includes
//
#include <iomanip>

// Project includes
//
#include "QuICC/Io/Variable/SphereTorPolNCoeffsPowerWriter.hpp"
#include "Environment/QuICCEnv.hpp"
#include "Types/Math.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"

namespace QuICC {

namespace Io {

namespace Variable {
   SphereTorPolNCoeffsPowerWriter::SphereTorPolNCoeffsPowerWriter(const std::string& prefix, const std::string& type, std::vector<std::shared_ptr<QuICC::DenseSM::Worland::RadialTorPolFunction>> pF)
      : ISphericalTorPolNScalarSpectrumWriter(prefix, type, std::vector<std::shared_ptr<QuICC::DenseSM::IGenericProfile>>(pF.begin(), pF.end()))
   {
   }

   SphereTorPolNCoeffsPowerWriter::~SphereTorPolNCoeffsPowerWriter() = default;

   void SphereTorPolNCoeffsPowerWriter::init()
   {
      // Sphere volume: 4/3*pi*r_o^3
      this->mVolume = (4.0 / 3.0) * Math::PI;

      this->mHasMOrdering = this->res().sim().ss().has(
         SpatialScheme::Feature::TransformSpectralOrdering123);

      ISphericalTorPolNScalarSpectrumWriter::init();
   }

}
}
}
